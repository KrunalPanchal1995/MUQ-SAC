/*  The functions in this file compute the broadening factors.
 *  The names of the functions start with "Fc" followed by
 *  the label of the corresponding reaction.
 */

#include"BroadeningFactors.h"

/*	Mechanism file: "Mech_li_et_al_2007.mech"	*/

BFFunction gBroadening[8] = { Fc13, Fc20, Fc26, Fc52, Fc53, Fc81, Fc82, Fc83 };

#ifndef MECHANISM
#define MECHANISM ""
#endif

void TReaction::CheckBroadeningFactors( const char *mechName )
{
	char	*name = new char[strlen( MECHANISM ) + 6];
	sprintf( name, "/%s.pre", MECHANISM );
	if ( strstr( mechName, name ) == NULL ) {
		for ( int i = 0; i < 8; ++i ) {
			gBroadening[i] = FcErr;
		}
	}
}

Double FcErr( Double /*T*/ )
{
	fprintf( stderr, "#error: wrong broadening factors (%s) linked to program\n", MECHANISM );
	exit( 2 );

	return 0;
}

Double Fc13( Double T )
{
#line 22 "Mech_li_et_al_2007.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc20( Double T )
{
#line 35 "Mech_li_et_al_2007.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc26( Double T )
{
#line 44 "Mech_li_et_al_2007.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc52( Double T )
{
#line 74 "Mech_li_et_al_2007.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc53( Double T )
{
#line 79 "Mech_li_et_al_2007.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc81( Double T )
{
#line 111 "Mech_li_et_al_2007.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc82( Double T )
{
#line 116 "Mech_li_et_al_2007.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc83( Double T )
{
#line 121 "Mech_li_et_al_2007.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

