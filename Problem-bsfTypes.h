/*==============================================================================
Project: Bulk Synchronous Farm (BSF)
Theme: BSF Skeleton
Module: Problem-bsfTypes.h (Predefined Problem-depended BSF Types)
Prefix: PT_bsf
Author: Leonid B. Sokolinsky
This source code is a part of BSF Skeleton
==============================================================================*/
#pragma once
#include "Problem-Types.h"		// Problem Types 

//=========================== BSF Types =========================
struct PT_bsf_parameter_T {			// Parameter for workers
	double approximation[PP_N];		// Current approximation
};

struct PT_bsf_mapElem_T {			// Element of map list
	int rowNo;						// Row number in matrix Alpha
};

struct PT_bsf_reduceElem_T {		// Element of reduce list	
	double g[PP_N];					// Coordinate
};

struct PT_bsf_reduceElem_T_1 {
	// optional filling
};

struct PT_bsf_reduceElem_T_2 {
	// optional filling
};

struct PT_bsf_reduceElem_T_3 {
	// optional filling
};