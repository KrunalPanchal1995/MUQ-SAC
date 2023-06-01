#define MECHANISM "Mech_li_et_al_2007"
#include "FlameMaster.h"

/*	Mechanism file: "Mech_li_et_al_2007.mech"	*/

typedef Double (*BFFunction)(Double T);

/* prototypes */
Double Fc13( Double T );
Double Fc20( Double T );
Double Fc26( Double T );
Double Fc52( Double T );
Double Fc53( Double T );
Double Fc81( Double T );
Double Fc82( Double T );
Double Fc83( Double T );
Double FcErr( Double T );


extern BFFunction gBroadening[8];
