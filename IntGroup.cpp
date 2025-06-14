#include <immintrin.h>
#include "IntGroup.h"
#include "Int.h"

using namespace std;

IntGroup::IntGroup(int size) {
  this->size = size;
  subp = (Int *)malloc(size * sizeof(Int));
}

IntGroup::~IntGroup() {
  free(subp);
}

void IntGroup::Set(Int *pts) {
  ints = pts;
}

// Compute modular inversion of the whole group
void IntGroup::ModInv() {

  Int newValue;
  Int inverse;

  subp[0].Set(&ints[0]);
  for (int i = 1; i < size; i++) {
    subp[i].ModMulK1(&subp[i - 1], &ints[i]);
  }

  // Do the inversion
  inverse.Set(&subp[size - 1]);
  inverse.ModInv();

  for (int i = size - 1; i > 0; i--) {
    newValue.ModMulK1(&subp[i - 1], &inverse);
    inverse.ModMulK1(&ints[i]);
    ints[i].Set(&newValue);
  }

  ints[0].Set(&inverse);

}

void IntGroup::ModInvK1order()
{
    if (size <= 0) return;               // пустой набор – делать нечего

    //-----------------------------------------------------------------
    // 1. Префикс-произведения:
    //    subp[i] =  a₀ · a₁ · … · aᵢ  (mod n)
    //-----------------------------------------------------------------
    subp[0].Set(&ints[0]);               // subp[0] = a₀
    for (int i = 1; i < size; ++i) {
        subp[i].Set(&subp[i - 1]);       // subp[i] = subp[i-1]
        subp[i].ModMulK1order(&ints[i]); //           · aᵢ   (mod n)
    }

    //-----------------------------------------------------------------
    // 2. Инверсия общего произведения
    //-----------------------------------------------------------------
    Int inverse(subp[size - 1]);         // inverse = (∏ aᵢ)
    inverse.ModInvK1();                  // inverse = (∏ aᵢ)⁻¹ (mod n)

    //-----------------------------------------------------------------
    // 3. Обратный проход: восстанавливаем каждое aᵢ⁻¹
    //-----------------------------------------------------------------
    for (int i = size - 1; i > 0; --i) {
        Int newValue(subp[i - 1]);       // subp[i-1] = a₀…aᵢ₋₁
        newValue.ModMulK1order(&inverse);// aᵢ⁻¹ = subp[i-1] · inverse

        inverse.ModMulK1order(&ints[i]); // inverse *= aᵢ  (готово для i-1)
        ints[i].Set(&newValue);          // записываем aᵢ⁻¹
    }
    ints[0].Set(&inverse);               // a₀⁻¹ осталось в «inverse»
}

