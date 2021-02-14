/*==============================================================================
Project: Bulk Synchronous Farm (BSF)
Theme: BSF Skeleton
Module: Problem-bsfCode.cpp (Problem-dependent Code)
Prefix: PI
Author: Leonid B. Sokolinsky
This source code is a part of BSF Skeleton
==============================================================================*/
#include "Problem-Data.h"			// Problem Types 
#include "Problem-Forwards.h"		// Problem Function Forwards
#include "Problem-bsfParameters.h"	// BSF-skeleton parameters
#include "BSF-SkeletonVariables.h"	// Skeleton Variables
using namespace std;

//----------------------- Predefined problem-dependent functions -----------------
void PC_bsf_Init(bool* success) { // success=false if initialization is unsuccessful
	PD_g = (double*)malloc(PP_N * sizeof(double));
	if (PD_g == NULL) {
		*success = false;
		return;
	}
}

void PC_bsf_SetListSize(int* listSize) {
	*listSize = PP_N;
};

void PC_bsf_MapF(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T* reduceElem, int* success // 1 - reduceElem was produced successfully (default); 0 - otherwise
){
	int mapIndex = BSF_sv_numberInSublist;

	if (mapIndex == 0)
		for (int i = 0; i < PP_N; i++)
			PD_g[i] = 0;

	for (int j = 0; j < PP_N; j++)
		PD_g[mapIndex + BSF_sv_addressOffset] += mapElem->row[j] * BSF_sv_parameter.approximation[j];

	if (mapIndex != BSF_sv_sublistLength - 1)
		*success = 0;
	else
		for (int j = 0; j < PP_N; j++)
			reduceElem->g[j] = PD_g[j];
};

void PC_bsf_MapF_1(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_1* reduceElem,
	int* success // 1 - reduceElem was produced successfully (default); 0 - otherwise
) {
	/* not used */
};

void PC_bsf_MapF_2(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_2* reduceElem,
	int* success // 1 - reduceElem was produced successfully (default); 0 - otherwise
) {
	/* not used */
};

void PC_bsf_MapF_3(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_3* reduceElem,
	int* success // 1 - reduceElem was produced successfully (default); 0 - otherwise
	) {
	// optional filling
};

void PC_bsf_ReduceF(PT_bsf_reduceElem_T* x, PT_bsf_reduceElem_T* y, PT_bsf_reduceElem_T* z) { // z = x + y
	for (int j = 0; j < PP_N; j++)
		z->g[j] = x->g[j] + y->g[j];
};

void PC_bsf_ReduceF_1(PT_bsf_reduceElem_T_1* x, PT_bsf_reduceElem_T_1* y, PT_bsf_reduceElem_T_1* z) {/* not used */ };
void PC_bsf_ReduceF_2(PT_bsf_reduceElem_T_2* x, PT_bsf_reduceElem_T_2* y, PT_bsf_reduceElem_T_2* z) {/* not used */ };
void PC_bsf_ReduceF_3(PT_bsf_reduceElem_T_3* x, PT_bsf_reduceElem_T_3* y, PT_bsf_reduceElem_T_3* z) {/* not used */ }

void PC_bsf_ProcessResults(
	PT_bsf_reduceElem_T* reduceResult,
	int reduceCounter, // Number of successfully produced Elrments of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
	for (int j = 0; j < PP_N; j++) {
		PD_prevApproximation[j] = parameter->approximation[j];
		parameter->approximation[j] = reduceResult->g[j] + beta(j);
	};

	/*debug*//*
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	cout << rank << ":=========>PC_bsf_ProcessResults: parameter->approximation:" << endl;
	for (int j = 0; j < PP_N; j++)
		cout << setw(8) << parameter->approximation[j];
	cout << endl;
	/*end debug*/

	if (ExitCondition(parameter))
		*exit = true;
	else
		*exit = false;
};

void PC_bsf_ProcessResults_1(
	PT_bsf_reduceElem_T_1* reduceResult,
	int reduceCounter, // Number of successfully produced Elrments of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
	// not used
};

void PC_bsf_ProcessResults_2(
	PT_bsf_reduceElem_T_2* reduceResult,
	int reduceCounter, // Number of successfully produced Elrments of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
	// not used
};

void PC_bsf_ProcessResults_3(
	PT_bsf_reduceElem_T_3* reduceResult,
	int reduceCounter, // Number of successfully produced Elrments of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
	// optional filling
};

void PC_bsf_JobDispatcher(
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* job,
	bool* exit
) {
	// Optional filling. Do not delete!
}

void PC_bsf_ParametersOutput(PT_bsf_parameter_T parameter) {
	cout << "=================================================== Jacobi M ====================================================" << endl;
	cout << "Number of Workers: " << BSF_sv_numOfWorkers << endl;
#ifdef PP_BSF_OMP
#ifdef PP_BSF_NUM_THREADS
	cout << "Number of Threads: " << PP_BSF_NUM_THREADS << endl;
#else
	cout << "Number of Threads: " << omp_get_num_procs() << endl;
#endif // PP_BSF_NUM_THREADS
#else
	cout << "OpenMP is turned off!" << endl;
#endif // PP_BSF_OMP

	cout << "Dimension: N = " << PP_N << endl;
	cout << "Eps_Square = " << PP_EPS << endl;
#ifdef PP_MATRIX_OUTPUT
	cout << "------- Matrix A & Column b -------" << endl;
	for (int i = 0; i < PP_N; i++) {
		for (int j = 0; j < PP_N; j++)
			cout << setw(7) << A(i, j);
		cout << setw(7) << b(i) << endl;
	};
	cout << endl;
	cout << "------- Matrix Alpha & Column Beta -------" << endl;
	for (int i = 0; i < PP_N; i++) {
		for (int j = 0; j < PP_N; j++)
			cout << setw(7) << Alpha(i, j);
		cout << setw(7) << beta(i) << endl;
	};
	cout << endl;
#endif // PP_MATRIX_OUTPUT
	cout << "Initial approximation: "; for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(7) << parameter.approximation[j]; cout << (PP_OUTPUT_LIMIT < PP_N ? "..." : "") << endl;
	cout << "-------------------------------------------" << endl;
};

void PC_bsf_CopyParameter(PT_bsf_parameter_T parameterIn, PT_bsf_parameter_T* parameterOutP) {
	for (int i = 0; i < PP_N; i++)
		parameterOutP->approximation[i] = parameterIn.approximation[i];
};

void PC_bsf_IterOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	cout << "------------------ " << BSF_sv_iterCounter << " ------------------" << endl;
	cout << "Approximation:\t\t"; for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(12) << parameter.approximation[j]; cout << (PP_OUTPUT_LIMIT < PP_N ? "..." : "") << endl;/**/

};

void PC_bsf_IterOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	cout << "------------------ " << BSF_sv_iterCounter << " ------------------" << endl;
	/* not used */
};

void PC_bsf_IterOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	cout << "------------------ " << BSF_sv_iterCounter << " ------------------" << endl;
	/* not used */
};

void PC_bsf_IterOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int jobCase) {
	cout << "------------------ " << BSF_sv_iterCounter << " ------------------" << endl;
	// optional filling

};

void PC_bsf_ProblemOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {// Output Function
	cout << "=============================================" << endl;
	cout << "Time: " << t << endl;
	cout << "Iterations: " << BSF_sv_iterCounter << endl;
	cout << "Solution: "; for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(12) << parameter.approximation[j]; cout << (PP_OUTPUT_LIMIT < PP_N ? "..." : "") << endl;/**/
};

void PC_bsf_ProblemOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {// Output Function
	// optional filling
};

void PC_bsf_ProblemOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {// Output Function
	// optional filling
};

void PC_bsf_SetInitParameter(PT_bsf_parameter_T* parameter) {
	for (int i = 0; i < PP_N; i++) // Generating coordinates of initial appriximation
		parameter->approximation[i] = beta(i);
};

void PC_bsf_ProblemOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double t) {// Output Function
	// optional filling
};

void PC_bsf_SetMapListElem(PT_bsf_mapElem_T* elem, int i) {
	for (int j = 0; j < PP_N; j++)
		elem->row[j] = Alpha(i, j);
}

//----------------------- Assigning Values to BSF-skeleton Variables (Do not modify!) -----------------------
void PC_bsfAssignAddressOffset(int value) { BSF_sv_addressOffset = value; };
void PC_bsfAssignIterCounter(int value) { BSF_sv_iterCounter = value; };
void PC_bsfAssignJobCase(int value) { BSF_sv_jobCase = value; };
void PC_bsfAssignMpiMaster(int value) { BSF_sv_mpiMaster = value; };
void PC_bsfAssignMpiRank(int value) { BSF_sv_mpiRank = value; };
void PC_bsfAssignNumberInSublist(int value) { BSF_sv_numberInSublist = value; };
void PC_bsfAssignNumOfWorkers(int value) { BSF_sv_numOfWorkers = value; };
void PC_bsfAssignParameter(PT_bsf_parameter_T parameter) { PC_bsf_CopyParameter(parameter, &BSF_sv_parameter); }
void PC_bsfAssignSublistLength(int value) { BSF_sv_sublistLength = value; };

//----------------------------- User functions -----------------------------
static double Norm(PT_vector_T x) {
	double sum = 0;
	for (int j = 0; j < PP_N; j++)
		sum += x[j] * x[j];
	return sqrt(sum);
};

static bool ExitCondition(PT_bsf_parameter_T* parameter) {
	PT_vector_T difference; // Difference between current and previous approximations

#ifdef PP_MAX_ITER_COUNT
	if (BSF_sv_iterCounter > PP_MAX_ITER_COUNT) {
		cout << "Acceptable maximum number of iterations is exceeded: PP_MAX_ITER_COUNT = " << PP_MAX_ITER_COUNT << endl;
		return true;
	};
#endif // PP_MAX_ITER_COUNT

	for (int j = 0; j < PP_N; j++)
		difference[j] = PD_prevApproximation[j] - parameter->approximation[j];
	if (Norm(difference) < PP_EPS)
		return true;
	else
		return false;
};

inline PT_float_T A(int i, int j) {
	if (j > i) return 0;
	if (j < i) return 1;
	if (j == 0) return 1;
	return 2 * i;
}

inline PT_float_T b(int i) {
	if (i == 0) return 1;
	return i + 2 * i;
}

inline PT_float_T Alpha(int i, int j) { // Reduced matrix Alpha
	if (j == i) return 0;
	return -A(i, j) / A(i, i);
}

inline PT_float_T beta(int i) { // Reduced colunm beta
	return b(i) / A(i, i);
}