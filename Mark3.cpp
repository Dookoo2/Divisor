// Mark3_updated.cpp  — дет-перебор “делителей” с батч-инверсией,
// + SIMD-bloom + DP-таблица + live-статистика
// + SHRINK-механизм: если решения нет за 10 минут,
//   secret ← secret − ⌊secret/20⌋; диапазон [1, secret].

#include <atomic>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include <omp.h>

#include "Int.h"
#include "Point.h"
#include "SECP256K1.h"
#include "IntGroup.h"
#include "simd_block_bloom.h"

//──────────────────── format helpers ───────────────────
static std::string humanBytes(size_t bytes){
    static const char* u[]={"B","KB","MB","GB","TB"};
    double v=double(bytes); int k=0;
    while(v>=1024.0&&k<4){ v/=1024.0; ++k; }
    std::ostringstream o;
    if(v<10)       o<<std::fixed<<std::setprecision(2);
    else if(v<100) o<<std::fixed<<std::setprecision(1);
    else           o<<std::fixed<<std::setprecision(0);
    o<<v<<u[k];
    return o.str();
}
static std::string formatRate(double cps){
    std::ostringstream o;
    if(cps>=1e6)      o<<std::fixed<<std::setprecision(1)<<cps/1e6<<"Mkey/s";
    else if(cps>=1e3) o<<std::fixed<<std::setprecision(1)<<cps/1e3<<"kkey/s";
    else              o<<std::fixed<<std::setprecision(0)<<cps<<"key/s";
    return o.str();
}

//──────────────────── SECP256K1 consts ───────────────────
static const char *N_HEX =
 "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141";
static const char *P_HEX =
 "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F";
static inline Int hexToInt(const char *h){
    Int x; x.SetBase16(const_cast<char*>(h)); return x;
}
static Int ORDER_N, FIELD_P;

//──────────────────── Scalar256 (для DP) ───────────────────
struct Scalar256{ uint64_t w[4]; };
static inline void intToScalar(const Int&src,Scalar256&dst){
    dst.w[0]=src.bits64[0]; dst.w[1]=src.bits64[1];
    dst.w[2]=src.bits64[2]; dst.w[3]=src.bits64[3];
}
static inline void scalarToInt(const Scalar256&s,Int&dst){
    dst.SetInt32(0);
    for(int i=3;i>=0;--i){ dst.ShiftL(64); Int t(s.w[i]); dst.Add(&t); }
}
static inline bool sameScalar(const Scalar256&a,const Scalar256&b){
    return std::memcmp(&a,&b,sizeof(a))==0;
}

//──────────────────── misc utils ───────────────────
static inline uint64_t splitmix64(uint64_t x){
    x+=0x9E3779B97F4A7C15ULL;
    x=(x^(x>>30))*0xBF58476D1CE4E5B9ULL;
    x=(x^(x>>27))*0x94D049BB133111EBULL;
    return x^(x>>31);
}
static inline uint64_t IntLow64(const Int&n){ return n.bits64[0]; }
static std::string toHex(const Int&v,bool pad=false){
    Int t; t.Set(&const_cast<Int&>(v));
    std::string s=t.GetBase16();
    for(char&c:s) c=std::tolower(c);
    if(pad&&s.size()<64) s.insert(0,64-s.size(),'0');
    return s;
}

//──────────────────── global objects ───────────────────
static Secp256K1 secp;
using fp_t = uint64_t;
static std::vector<fp_t>      fp_tbl;
static std::vector<Scalar256> idx_tbl;
static std::unique_ptr<std::atomic<uint8_t>[]> used_tbl;
static uint32_t dp_cap = 0;
static simd_bloom::SimdBlockFilterFixed<> *bloom = nullptr;

//──────────────────── Bloom + DP ───────────────────
static inline fp_t make_fp(const Point&P){
    return splitmix64(IntLow64(P.x)^uint64_t(!P.y.IsEven()));
}
static bool dp_insert_unique(fp_t fp,const Int&idx){
    Scalar256 key; intToScalar(idx,key);
    uint32_t h=uint32_t(fp)%dp_cap;
    for(;;){
        uint8_t st=used_tbl[h].load(std::memory_order_acquire);
        if(st==2){
            if(fp_tbl[h]==fp&&sameScalar(idx_tbl[h],key)) return false;
        }else if(st==0){
            uint8_t z=0;
            if(used_tbl[h].compare_exchange_strong(z,1,std::memory_order_acq_rel)){
                fp_tbl[h]=fp; idx_tbl[h]=key;
                used_tbl[h].store(2,std::memory_order_release);
                return true;
            }
        }
        if(++h==dp_cap) h=0;
    }
}
static bool dp_find(fp_t fp,Int&out){
    uint32_t h=uint32_t(fp)%dp_cap;
    while(used_tbl[h].load(std::memory_order_acquire)==2){
        if(fp_tbl[h]==fp){ scalarToInt(idx_tbl[h],out); return true; }
        if(++h==dp_cap) h=0;
    }
    return false;
}

//──────────────────── CLI args ───────────────────
struct Args{ std::string keyHex,rangeStr; uint64_t q=0; std::string order="forward"; } args;
static void usage(const char*p){
    std::cerr<<"Usage: "<<p<<" -k <pubkey_hex> -r <start_hex>:<end_hex>"
             <<" [-q <N>] [--order forward|random]\n";
}
static bool parseRange(const std::string&s,Int&a,Int&b){
    auto p=s.find(':'); if(p==std::string::npos) return false;
    a.SetBase16(const_cast<char*>(s.substr(0,p).c_str()));
    b.SetBase16(const_cast<char*>(s.substr(p+1).c_str()));
    return !a.IsGreater(&b);
}
static bool parseArgs(int ac,char**av){
    if(ac<5){ usage(av[0]); return false; }
    for(int i=1;i<ac;++i){
        std::string a=av[i];
        if(a=="-k") args.keyHex=av[++i];
        else if(a=="-r") args.rangeStr=av[++i];
        else if(a=="-q") args.q=std::strtoull(av[++i],nullptr,10);
        else if(a=="--order") args.order=av[++i];
        else{ std::cerr<<"Unknown arg "<<a<<"\n"; return false; }
    }
    if(args.keyHex.empty()||args.rangeStr.empty()){ usage(av[0]); return false; }
    if(args.order!="forward"&&args.order!="random"){
        std::cerr<<"--order must be forward|random\n"; return false;
    }
    return true;
}

//──────────────────── randomInRange ───────────────────
static Int randomInRange(const Int&start,const Int&range,std::mt19937_64&rng){
    Int r; for(int i=0;i<4;++i) r.bits64[i]=rng();
    r.Mod(const_cast<Int*>(&range));
    Int t(start); r.Add(&t); return r;
}

//──────────────────── SHRINK helpers ───────────────────
static std::mutex shrink_mtx;
static int        shrinkCount = 0;
static Int        offset;
static std::atomic<bool> time_up(false);

// Новый mutex и текущий счётчик делителей
static std::mutex cur_mtx;
static Int        global_cur;

static inline Point negatePoint(const Point&P){
    Point R=P;
    Int yn(FIELD_P), yCopy(P.y);
    yn.Sub(&yCopy); yn.Mod(&FIELD_P);
    R.y = yn;
    return R;
}
static inline void pointSub(Point&P,const Point&Q){
    Point neg = negatePoint(Q);
    P = secp.Add(P, neg);
}
static Int computeShrinkStep(const Int&curRange){
    Int twenty; twenty.SetInt32(10);
    Int step(curRange), rem;
    step.Div(&twenty,&rem);
    if(step.IsZero()) step.SetInt32(1);
    return step;
}

//──────────────────────── MAIN ─────────────────────────
int main(int argc,char**argv){
    if(!parseArgs(argc,argv)) return 1;

    ORDER_N = hexToInt(N_HEX);
    FIELD_P = hexToInt(P_HEX);
    secp.Init();
    offset.SetInt32(0);

    //── public key P
    Point P;
    if(args.keyHex.size()==66&&(args.keyHex[1]=='2'||args.keyHex[1]=='3')){
        bool even = (args.keyHex[1]=='2');
        Int x; x.SetBase16(const_cast<char*>(args.keyHex.substr(2).c_str()));
        Int y = secp.GetY(x, even);
        P.x=x; P.y=y; P.z.SetInt32(1);
    } else if(args.keyHex.size()==130&&args.keyHex[1]=='4'){
        Int x; x.SetBase16(const_cast<char*>(args.keyHex.substr(2,64).c_str()));
        Int y; y.SetBase16(const_cast<char*>(args.keyHex.substr(66,64).c_str()));
        P.x=x; P.y=y; P.z.SetInt32(1);
    } else { std::cerr<<"Bad key\n"; return 1; }

    //── начальный диапазон
    Int tmpStart, tmpEnd;
    if(!parseRange(args.rangeStr, tmpStart, tmpEnd)){ std::cerr<<"Bad range\n"; return 1; }
    Int one; one.SetInt32(1);
    Int start = one;
    Int end   = tmpEnd;
    Int curRange = end;

    //── Phase-0: allocate DP/Bloom
    if(args.q==0) args.q=1;
    const uint64_t traps=args.q;
    const double LOAD=0.75, BLOOM_F=2.0;
    dp_cap = uint32_t(std::ceil(traps/LOAD));

    fp_tbl.assign(dp_cap,0);
    idx_tbl.assign(dp_cap,Scalar256{0,0,0,0});
    used_tbl.reset(new std::atomic<uint8_t>[dp_cap]);
    for(uint32_t i=0;i<dp_cap;++i) used_tbl[i].store(0);
    bloom = new simd_bloom::SimdBlockFilterFixed<>(size_t(traps*BLOOM_F));

    size_t slotBytes  = sizeof(fp_t)+sizeof(Scalar256)+1;
    size_t dpBytes    = size_t(dp_cap)*slotBytes;
    size_t bloomBytes = size_t(traps*BLOOM_F);
    std::cout<<"\n=========== Phase-0: RAM summary ===========\n"
             <<"DP table : "<<humanBytes(dpBytes)<<"  ("<<traps<<" / "<<dp_cap
             <<" slots, load "<<std::fixed<<std::setprecision(2)
             <<double(traps)/dp_cap*100<<"%)\n"
             <<"Bloom    : "<<humanBytes(bloomBytes)<<"\n"
             <<"Total    : "<<humanBytes(dpBytes+bloomBytes)<<"\n\n";

    //── fill DP/Bloom
    Int mul; mul.SetInt32(1);
    for(uint64_t m=1; m<=traps; ++m){
        Point Q = secp.ComputePublicKey(&mul);
        fp_t fp = make_fp(Q);
        dp_insert_unique(fp, mul);
        bloom->Add(uint32_t(fp));
        mul.Add(&one);
    }

    //── warm-up fixed-base
    { Int d; d.SetInt32(1); secp.ComputeFixedPointMul(&d,P); }

    //── live-monitor
    std::atomic<bool> found(false);
    std::atomic<long long> globalChecks(0);
    auto t0 = std::chrono::high_resolution_clock::now();
    std::thread monitor([&]{
        while(!found.load()){
            std::this_thread::sleep_for(std::chrono::seconds(2));
            auto now = std::chrono::high_resolution_clock::now();
            double dt = std::chrono::duration<double>(now-t0).count();
            long long chk = globalChecks.load();
            double cps = chk/dt;
            std::lock_guard<std::mutex> lk(shrink_mtx);
            std::cout<<"\033[4A"
                     <<"Full time: "<<std::fixed<<std::setprecision(2)<<dt<<" s\n"
                     <<"Speed    : "<<formatRate(cps)<<"\n"
                     <<"Shrink   : "<<shrinkCount<<"\n"
                     <<"Offset   : 0x"<<toHex(offset,true)<<"\n";
        }
    });
    std::cout<<"============= Phase-1: Working =============\n"
             <<"Full time: 0.00 s\n"
             <<"Speed    : "<<formatRate(0.0)<<"\n"
             <<"Shrink   : 0\n"
             <<"Offset   : 0x0\n";

    //── result vars
    Int foundDiv, foundMul, privKey;

    //──────────────── outer SHRINK loop ────────────────
    while(curRange.IsGreater(&one) && !found.load()){
        // сбросить счётчик текущего делителя на 1
        {
            std::lock_guard<std::mutex> lk(cur_mtx);
            global_cur.SetInt32(1);
        }

        auto iterStart = std::chrono::high_resolution_clock::now();
        auto deadline  = iterStart + std::chrono::seconds(60);
        time_up.store(false);

        //──────── main search (OMP) ────────
        #pragma omp parallel
        {
            constexpr int BATCH = 256;
            std::vector<Int> orig(BATCH), buf(BATCH);
            IntGroup grp(BATCH); grp.Set(buf.data());

            while(!found.load() && !time_up.load()){
                // зарезервировать следующий батч [batch_start … batch_end]
                Int batch_start, batch_end;
                {
                    std::lock_guard<std::mutex> lk(cur_mtx);
                    if(global_cur.IsGreater(&curRange)) break;
                    batch_start = global_cur;
                    batch_end   = global_cur;
                    // вычисляем batch_end = min(global_cur + (BATCH-1), curRange)
                    for(int i=1; i<BATCH; ++i){
                        Int tmp = batch_end; tmp.ModAdd(&tmp, &one);
                        if(tmp.IsGreater(&curRange)) break;
                        batch_end = tmp;
                    }
                    // сдвинуть global_cur сразу за batch_end
                    global_cur = batch_end; global_cur.ModAdd(&global_cur, &one);
                }

                // заполнить orig/ buf последовательными значениями
                int cnt = 0;
                Int cur = batch_start;
                while(cnt < BATCH && !cur.IsGreater(&batch_end)){
                    orig[cnt] = cur;
                    buf[cnt]  = cur;
                    ++cnt;
                    cur.ModAdd(&cur, &one);
                }
                // если недобрали до BATCH, дублируем последний
                for(int i=cnt; i<BATCH; ++i){
                    orig[i] = orig[cnt-1];
                    buf[i]  = buf[cnt-1];
                }

                // проверка батча
                grp.ModInvK1order();
                for(int i=0; i<BATCH && !time_up.load(); ++i){
                    globalChecks.fetch_add(1);
                    Point Q = secp.ComputeFixedPointMul(&buf[i], P);
                    fp_t fp = make_fp(Q);
                    if(!bloom->Find(uint32_t(fp))) continue;
                    Int mulCand;
                    if(!dp_find(fp, mulCand)) continue;

                    Int priv(mulCand);
                    priv.ModMulK1(&priv, &orig[i]);
                    priv.Mod(&ORDER_N);

                    if(!found.exchange(true)){
                        #pragma omp critical
                        {
                            foundDiv = orig[i];
                            foundMul = mulCand;
                            privKey  = priv;
                        }
                    }
                    break;
                }

                if(std::chrono::high_resolution_clock::now() >= deadline)
                    time_up.store(true);
            }
        } // конец #pragma omp parallel

        if(found.load()) break;

        //──────── SHRINK ────────
        Int step = computeShrinkStep(curRange);
        {
            std::lock_guard<std::mutex> lk(shrink_mtx);
            ++shrinkCount;
            offset.Add(&step);
        }
        // P = P − step·G
        Point deltaG = secp.ComputePublicKey(&step);
        pointSub(P, deltaG);
        P.Reduce();

        // пересбор фиксированной базы
        secp.ClearFixedBase();
        { Int oneTmp; oneTmp.SetInt32(1); secp.ComputeFixedPointMul(&oneTmp, P); }

        // уменьшить диапазон
        curRange.Sub(&step);
        end = curRange;
    }

    //── завершение мониторинга
    found.store(true);
    monitor.join();

    //── Phase-2
    std::cout<<"\n============== Phase-2: RESULT =============\n";
    if(found){
        Int finalPriv(privKey);
        finalPriv.Add(&offset);
        finalPriv.Mod(&ORDER_N);
        std::cout<<"Divisor  : 0x"<<toHex(foundDiv,true)<<"\n"
                 <<"Multiple : 0x"<<toHex(foundMul,true)<<"\n"
                 <<"Offset   : 0x"<<toHex(offset,true)<<"\n"
                 <<"Private  : 0x"<<toHex(finalPriv,true)<<"\n";
    } else {
        std::cout<<"Not found (range exhausted)\n";
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout<<"Full time: "
             <<std::chrono::duration<double>(t1-t0).count()<<" s\n";

    delete bloom;
    return 0;
}
