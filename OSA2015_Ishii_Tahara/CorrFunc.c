/*	CorrFunc.c
*/

#include "XOPStandardHeaders.h"			// Include ANSI headers, Mac headers, IgorXOP.h, XOP.h and XOPSupport.h

// Prototypes
HOST_IMPORT void main(IORecHandle ioRecHandle);

// Custom error codes
#define REQUIRES_IGOR_400 1 + FIRST_XOP_ERR
//#define NON_EXISTENT_WAVE 2 + FIRST_XOP_ERR
//#define REQUIRES_SP_OR_DP_WAVE 3 + FIRST_XOP_ERR

#include "XOPStructureAlignmentTwoByte.h"	// All structures passed to Igor are two-byte aligned.

typedef struct calc2DCorrHistogramParams {
	double dtau;		// window size
	double tau;		// delay time
	double n2;		// length of the table 2
	double n1;		// length of the table 1
	waveHndl microtime2; // microtime table 2
	waveHndl microtime1; // microtime table 1
	waveHndl timetable2;	// macrotime table 2
	waveHndl timetable1;	// macrotime table 1
	double result;
} calc2DCorrHistogramParams, *calc2DCorrHistogramParamsPtr;

typedef struct calcPhotonAssocDecayParams {
	double dtau;		// window size
	double tau;		// delay time
	double n2;		// length of the table 2
	double n1;		// length of the table 1
	waveHndl microtime1; // microtime table 1
	waveHndl timetable2;	// macrotime table 2
	waveHndl timetable1;	// macrotime table 1
	double result;
} calcPhotonAssocDecayParams, *calcPhotonAssocDecayParamsPtr;

#include "XOPStructureAlignmentReset.h"

static int
calc2DCorrHistogram(calc2DCorrHistogramParamsPtr p)
{
	unsigned long n1,n2;
	unsigned long tau,dtau;
	
	waveHndl wavH;
	unsigned long *mp,*mp_offset;
	long dims[MAX_DIMENSIONS+1];
	int err;

	unsigned long i1,j1,j2;
	unsigned char idelay,jdelay;
	unsigned long t1,t2,t3,t4;
	unsigned long *tp1,*tp2;
	unsigned char *dp1,*dp2,*dp2_offset;
	int hState1,hState2,hState3,hState4,hState5;

	n1 = (unsigned long)p->n1;
	n2 = (unsigned long)p->n2;
	tau = (unsigned long)p->tau;
	dtau = (unsigned long)p->dtau;

	if(!((n1==WavePoints(p->timetable1))&&(n2==WavePoints(p->timetable2))
		&&(n1==WavePoints(p->microtime1))&&(n2==WavePoints(p->microtime2)))) {
			SetNaN64(&p->result);
			return 0;
		}

	if(!((WaveType(p->timetable1) & (NT_I32|NT_UNSIGNED))&&(WaveType(p->timetable2) & (NT_I32|NT_UNSIGNED))
		&&(WaveType(p->microtime1) & (NT_I8|NT_UNSIGNED))&&(WaveType(p->microtime2) & (NT_I8|NT_UNSIGNED)))) {
		SetNaN64(&p->result);
		return 0;
	}
	
	hState1=MoveLockHandle(p->timetable1);
	hState2=MoveLockHandle(p->timetable2);
	hState3=MoveLockHandle(p->microtime1);
	hState4=MoveLockHandle(p->microtime2);

	tp1 = (unsigned long*)WaveData(p->timetable1);
	tp2 = (unsigned long*)WaveData(p->timetable2);
	dp1 = (unsigned char*)WaveData(p->microtime1);
	dp2 = (unsigned char*)WaveData(p->microtime2);

	MemClear(dims, sizeof(dims));
	dims[ROWS] = 256;
	dims[COLUMNS] = 256;
	if (err = MDMakeWave(&wavH, "M_CorrHist", NULL, dims, NT_I32|NT_UNSIGNED, 1))
		return err;
	hState5=MoveLockHandle(wavH);
	mp = (unsigned long*)WaveData(wavH);
	MemClear(mp, 256*256*sizeof(unsigned long));
	
	j1 = 0;
	j2 = 0;
	for(i1=0;i1<n1;i1++) {
		t1 = tp1[i1];
		t3 = t1 + tau;
		t4 = t3 + dtau;
		if(t3 > tp2[n2-1])
			break;
		idelay = dp1[i1];
		mp_offset = mp + (idelay<<8);
		while((j1 < n2) && (tp2[j1] < t3))
			j1++;
		while((j2 < n2) && (tp2[j2] < t4))
			j2++;
		for(dp2_offset=dp2+j1;dp2_offset<dp2+j2;dp2_offset++)
			(*(mp_offset+*dp2_offset))++;
	}

	if(i1 == 0) {
		SetNaN64(&p->result);
		return 0;
	}

	HSetState(wavH, hState5);
	HSetState(p->microtime2, hState4);
	HSetState(p->microtime1, hState3);
	HSetState(p->timetable2, hState2);
	HSetState(p->timetable1, hState1);

	return 0;

}

static int
calcPhotonAssocDecay(calcPhotonAssocDecayParamsPtr p)
{
	unsigned long n1,n2;
	unsigned long tau,dtau;
	
	waveHndl wavH;
	unsigned long *mp,*mp_offset;
	long dims[MAX_DIMENSIONS+1];
	int err;

	unsigned long i1,j1,j2;
	unsigned char idelay,jdelay;
	unsigned long t1,t2,t3,t4;
	unsigned long *tp1,*tp2;
	unsigned char *dp1,*dp2,*dp2_offset;
	int hState1,hState2,hState3,hState4;

	n1 = (unsigned long)p->n1;
	n2 = (unsigned long)p->n2;
	tau = (unsigned long)p->tau;
	dtau = (unsigned long)p->dtau;

	if(!((n1==WavePoints(p->timetable1))&&(n2==WavePoints(p->timetable2))
		&&(n1==WavePoints(p->microtime1)))) {
			SetNaN64(&p->result);
			return 0;
		}

	if(!((WaveType(p->timetable1) & (NT_I32|NT_UNSIGNED))&&(WaveType(p->timetable2) & (NT_I32|NT_UNSIGNED))
		&&(WaveType(p->microtime1) & NT_I8|NT_UNSIGNED))) {
		SetNaN64(&p->result);
		return 0;
	}
	
	hState1=MoveLockHandle(p->timetable1);
	hState2=MoveLockHandle(p->timetable2);
	hState3=MoveLockHandle(p->microtime1);

	tp1 = (unsigned long*)WaveData(p->timetable1);
	tp2 = (unsigned long*)WaveData(p->timetable2);
	dp1 = (unsigned char*)WaveData(p->microtime1);

	MemClear(dims, sizeof(dims));
	dims[ROWS] = 256;
	if (err = MDMakeWave(&wavH, "W_PhotonAssocDecay", NULL, dims, NT_I32|NT_UNSIGNED, 1))
		return err;
	hState4=MoveLockHandle(wavH);
	mp = (unsigned long*)WaveData(wavH);
	MemClear(mp, 256*sizeof(unsigned long));
	
	j1 = 0;
	j2 = 0;
	for(i1=0;i1<n1;i1++) {
		t1 = tp1[i1];
		t3 = t1 + tau;
		t4 = t3 + dtau;
		idelay = dp1[i1];
		mp_offset = mp + idelay;
		while((j1 < n2) && (tp2[j1] < t3))
			j1++;
		while((j2 < n2) && (tp2[j2] < t4))
			j2++;
		*mp_offset += j2 - j1;
	}

	if(i1 == 0) {
		SetNaN64(&p->result);
		return 0;
	}

	HSetState(wavH, hState4);
	HSetState(p->microtime1, hState3);
	HSetState(p->timetable2, hState2);
	HSetState(p->timetable1, hState1);

	return 0;

}

static long
RegisterFunction()
{
	int funcIndex;

	funcIndex = GetXOPItem(0);			// Which function invoked ?
	switch (funcIndex) {
		case 0:							// y = SimpleFit(w,x) (curve fitting function).
			return((long)calc2DCorrHistogram);	// This function is called using the direct method.
			break;
		case 1:
			return((long)calcPhotonAssocDecay);
			break;
	}
	return NIL;
}

/*	XOPEntry()

	This is the entry point from the host application to the XOP for all
	messages after the INIT message.
*/
static void
XOPEntry(void)
{	
	long result = 0;

	switch (GetXOPMessage()) {
		case FUNCADDRS:
			result = RegisterFunction();	// This tells Igor the address of our function.
			break;
	}
	SetXOPResult(result);
}

/*	main(ioRecHandle)

	This is the initial entry point at which the host application calls XOP.
	The message sent by the host must be INIT.
	main() does any necessary initialization and then sets the XOPEntry field of the
	ioRecHandle to the address to be called for future messages.
*/
HOST_IMPORT void
main(IORecHandle ioRecHandle)
{	
	XOPInit(ioRecHandle);							// Do standard XOP initialization.
	SetXOPEntry(XOPEntry);							// Set entry point for future calls.
	
	if (igorVersion < 400)
		SetXOPResult(REQUIRES_IGOR_400);
	else
		SetXOPResult(0L);
}
