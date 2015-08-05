#ifndef __clasa__
#define __clasa__


#include <iostream>
#include "TObject.h"


class clasa : public TObject {
    int i;
    public:
    void seti(int);
    int geti();
    clasa(){i=0;};
    ~clasa();
    ClassDef(clasa, 1);
};

#endif
