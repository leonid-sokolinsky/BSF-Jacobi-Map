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
	PT_vector_T approximation;		// Current approximation
};

struct PT_bsf_mapElem_T {		// Element of map list
	PT_vector_T row;			// Row of reduced matrix		
};

struct PT_bsf_reduceElem_T {		// Element of reduce list	
	PT_vector_T g;					// Coordinate
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