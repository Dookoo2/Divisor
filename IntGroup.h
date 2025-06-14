#ifndef INTGROUPH
#define INTGROUPH

#include "Int.h"
#include <vector>

class IntGroup {

public:

	IntGroup(int size);
	~IntGroup();
    void Set(Int* pts);                // массив входных чисел
    void ModInv();                     // старая:   по модулю P  (кривой)
    void ModInvK1order();              // новая:   по порядку n  (группа)

private:

    int  size {0};
    Int* ints {nullptr};               // массив входных чисел
    Int* subp {nullptr};               // промежуточные префикс-произведения

};

#endif // INTGROUPCPUH
