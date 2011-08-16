// cnF2freq, (c) Carl Nettelblad, Department of Information Technology, Uppsala University 2008
// Public release 0.5
//
// carl.nettelblad@it.uu.se
//
// This code is allowed to be freely used for any commercial or research purpose. If the code is
// used or integrated into another project largely unchanged, attribution to the original author
// is appreciated. No warranties are given.


#define NDEBUG
// These defines fixed an error in one particular site installation of the Portland compiler.
#define _STLP_EXPOSE_GLOBALS_IMPLEMENTATION 1
#define _REENTRANT 1
#define _SECURE_SCL 0
// For MSCVC
// Note: MSVC OpenMP support is insufficient in current release. Disable OpenMP for compilation
// in MSVC.
//
// Recent releases of g++ on MacOS X and Linux, as well as the Intel C++ compiler on
// Linux (x86 and x64) and Windows have been tested. The Portland compiler collection works only
// with some optimization settings, some race conditions in STL, although some
// workarounds are used.
#define _CRT_SECURE_NO_WARNINGS

#include <vector>

#include <stdio.h>

float templgeno[8] = {-1, -0.5,
	0,  0.5,
	0,  0.5,
	1, -0.5};



// _MSC_VER is here to be interpreted as any compiler providing TR1 C++ headers
#ifdef _MSC_VER
#include <array>
#else
// Boost also provides an array implementation, that is largely compatible
#include <boost/array.hpp>
#endif

#include "cnf2freqmini.h"
#include <boost/static_assert.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>




#include <assert.h>
#include <stdlib.h>
#include <set>
#include <algorithm>
#include <math.h>
#include <map>
#include <float.h>
#include <string.h>


using namespace std; // use functions that are part of the standard library
#ifdef _MSC_VER
using namespace tr1;
#else
using namespace boost;
#endif

// Define needed due to namespace pollution
#define none cnF2freqNONE


#ifndef _MSC_VER
#define _isnan isnan
#define _finite isfinite
#endif

// Data structures for markers/chromosomes
std::vector<double> markerposes;
std::vector<double> actrec[2];
std::vector<unsigned int> chromstarts;
std::vector<bool> sexmarker;



// Special casting considering negative values only to be true
struct negtrue
{
	const int value;
	negtrue(bool neg, int value) : value(neg ? -1 : (value >= 0 ? value : 0))
	{	
	}

	operator bool() const
	{
		return value < 0;
	}
};


#ifndef MCIBD
// F2 with no haplotyping
const int NUMGEN = 2;
const int TYPEBITS = (1 << NUMGEN) - 2;
const int TYPESEXES[TYPEBITS] = {0, 1};
const int NUMTYPES = 1 << TYPEBITS;
const double EVENGEN = 1.0 / NUMTYPES;
const float MINFACTOR = -1e15;
const unsigned int NUMFLAG2GEN = 1;
const unsigned int HALFNUMPATHS = 1;
const unsigned int NUMPATHS = 2;
const unsigned int NUMSHIFTGEN = 0;
const unsigned int HALFNUMSHIFTS = 1;
const unsigned int NUMSHIFTS = 1;
const bool HAPLOTYPING = false;
#else
const int NUMGEN = 3;
const int TYPEBITS = (1 << NUMGEN) - 2;
const int TYPESEXES[TYPEBITS] = {0, 0, 1, 1, 0, 1};
const int NUMTYPES = 1 << TYPEBITS;
const double EVENGEN = 1.0 / NUMTYPES;
const float MINFACTOR = -1e15;
const unsigned int NUMFLAG2GEN = NUMGEN;
const unsigned int HALFNUMPATHS = 1 << (TYPEBITS / 2);
const unsigned int NUMPATHS = NUMTYPES << 1;
const unsigned int NUMSHIFTGEN = 1;
const unsigned int HALFNUMSHIFTS = 1;
const unsigned int NUMSHIFTS = 1; 
const bool HAPLOTYPING = false;
#endif


int homogamsex = -1;


const int HALFNUMTYPES = 1 << (TYPEBITS / 2);

// Infer corrections for impossible genotypes according to the pedigree
// Also infer values for missing markers from existing information in offspring
// When enabled, results similar to QTL Express without data reduction
// ccoeff does not provide correction inference, so exact result reproduction
// is achieved when this flag is disabled.
bool CORRECTIONINFERENCE = false;

bool early = false;


typedef vector<array<double, 5 > > MWTYPE;


MWTYPE markerweight;

double discstep = 0.1;
double baserec[2];
int sexc = 1;


// Is a specific marker value a admissible as a match to marker b
// The logic we use currently represents "0" as unknown, anything else as a known
// marker value.
// NOTE, VALUE CHANGED FOR ZEROS!
template<bool zeropropagate> bool markermiss(MarkerVal& a, const MarkerVal b)
{
	// This is the real logic; we do not bind anything at all when zeropropagate is true
	if (zeropropagate) return false;

	if (a == UnknownMarkerVal)
	{
		if (!zeropropagate) a = b;
		return false;
	}
	if (b == UnknownMarkerVal && a != sexmarkerval) return false;

	return a != b;
}


// Move a tree-based binary flag up a generation. The structure of bit flags might look like
// (1)(3)(3), where each group of (3) is in itself (1)(2), where (2) is of course (1)(1).
int upflagit(int flag, int parnum, int genwidth)
{

	if (flag < 0) return flag;
	flag >>= parnum * (genwidth - 1);
	flag &= ((1 << (genwidth - 1)) - 1);

	return flag;
}

struct individ;

int generation = 1;
int shiftflagmode;


// Convenient helpers to represent arrays and matrices connected to our number of states
template<class T2> class PerStateArray
{
public:
	typedef array<T2, NUMTYPES> T;
};

template<class T2> class StateToStateMatrix
{
public:
	typedef typename PerStateArray<
		typename PerStateArray<T2>::T >
		::T T;
};


// done, factors and cacheprobs all keep track of the same data
// done indicates that a specific index (in the binary tree of blocks of multi-step transitions) is done
// with a "generation id" that's semi-unique, meaning no active clearing of the data structure is performed
vector<int> done[NUMSHIFTS];
// factors contain the mantissas of the extended floating-point representation
vector<PerStateArray<float>::T > factors[NUMSHIFTS];
// cacheprobs contain actual transitions from every possible state to every possible other state
vector<StateToStateMatrix<float>::T > cacheprobs[NUMSHIFTS];
vector<individ*> reltree;


// We put all thread-local structures in our own separate struct. This is because many compilers implementing OpenMP
// use a relatively expensive call for determining thread-local global data, which is not cached over function calls.
// The simple pointer-arithmetic for a lookup within a struct is cheaper.
struct threadblock
{
	int* const generation;
	int* const shiftflagmode;


	vector<int>* const done;
	vector<PerStateArray<float>::T >* const factors;
	vector<StateToStateMatrix<float>::T >* const cacheprobs;	

	threadblock() : generation(&::generation), shiftflagmode(&::shiftflagmode),
		done(::done), factors(::factors), cacheprobs(::cacheprobs)
		
	{
	};
};

// Turners are used as mix-ins when probabilities are filtered at the "fixed" marker
// (the one we are actually probing).
// This functionality is used to invert the state in different manner, to test for the
// possibility that the complete haplotype assignment from an arbitrary point until the end
// has been mixed up. This can easily happen as the assignment of haplotype numbers is
// basically arbitrary, so the only thing keeping it in place is the linkage to the original
// defining locus.
class noneturner
{
public:
	void operator() (const PerStateArray<double>::T& probs) const
	{
	};

	bool canquickend() const
	{
		return true;
	}
} none;


// A struct containing some auxiliary arguments to the trackpossible family of functions.
// These parameters are not known at compile-time, hence not template arguments, but they do not
// change over the series of recursive calls.
const struct trackpossibleparams
{
	const float updateval;
	int* const gstr;

	trackpossibleparams() : updateval(0.0f), gstr(0)
	{
	}

	trackpossibleparams(float updateval, int* gstr) : updateval(updateval), gstr(gstr)
	{
	}
} tpdefault;

struct classicstop
{
	int lockpos;
	int genotype;

	classicstop(int lockpos, int genotype) : lockpos(lockpos), genotype(genotype)
	{}

	operator int() const
	{
		return lockpos;
	}

	const int getgenotype(int startmark) const
	{
		return genotype;
	}

	negtrue okstep(const int startmark, const int endmark) const
	{
		bool found = (lockpos <= -1000 - startmark && lockpos > -1000 - endmark) ||
			(lockpos >= markerposes[startmark] && lockpos <= markerposes[endmark]);

		return negtrue(!found, genotype);
	}

	bool fixtofind(int& genotype, double& startpos, double& endpos, int j) const
	{
		bool tofind = (lockpos >= startpos && lockpos <= endpos);
		if (tofind)
		{
			endpos = lockpos;
		}


		if (lockpos == -1000 - j + 1)
		{
			tofind = true;
			endpos = startpos;
		}

		if (tofind)
		{
			genotype = this->genotype;
		}

		return tofind;
	}

};


struct smnonecommon
{
	const int getgenotype(int startmark) const
	{
		return -1;
	}
};

struct nonestop : smnonecommon
{
	operator int() const
	{
		assert(false);

		return -1;
	}

	negtrue okstep(const int startmark, const int endmark) const
	{
		return negtrue(true, 0);
	}

	bool fixtofind(int& genotype, double& startpos, double& endpos, int j) const
	{
		return false;
	}
} NONESTOP;


template<class G> const bool canquickend(int startmark, const G &stopdata)
{
	return false;
}

/*template<> const bool canquickend<smnonecommon>(int startmark, const smnonecommon &stopdata)
{
return false;
}*/

template<> const bool canquickend<classicstop>(int startmark, const classicstop &stopdata)
{
	return stopdata.lockpos <= -1000 && stopdata.genotype >= 0;
}



// Should only be called for markers where locking is really done
template<class G> bool firstlockmatch(int lockpos, const G& stopdata)
{
	return stopdata.lockpos == lockpos;
}

template<> bool firstlockmatch<nonestop>(const int lockpos, const nonestop& stopdata)
{
	return false;
}


template<class G> double getactrec(const G& stopdata, double startpos, double endpos, int k, int j)
{
	return actrec[k][j + 1];
}



// A structure containing most of the information on an individual
struct individ : miniindividual
{	
	// Temporary storage of all possible marker values, used in fixparents.
	vector<set<MarkerVal> > markervals;

	//
	vector<individ*> kids;

	// Number of children.
	int children;
	individ** parsalias;
	individ**& pars;

	individ() : parsalias((individ**) &(minipars)), pars(parsalias)
	{
		pars[0] = 0;
		pars[1] = 0;
		strain = 0;
		sex = false;
	}



	// A wrapper for calling trackpossible. Additional logic introduced to do lookup in the "impossible" tables of branches
	// not allowed at the current point. This will generally mean most branches and reduces the number of states visited considerably
	// in typical cases.
	//
	// The double overload means that the class can be treated as a conventional function call, returning a double, when the pre-lookup
	// is not needed. In the "real" trackpossible method, one instance is called in that way, while the other one is pre-looked up.
	template<bool update, bool zeropropagate> struct recursetrackpossible
	{
		individ* mother;
		int upflagr;
		int upflag2r;
		int upshiftr;
		const trackpossibleparams& extparams;
		const int genwidth;
		const MarkerVal markerval;
		const int marker;
		const threadblock& tb;
		int firstpar;
		int* impossibleref;
		int impossibleval;
		bool prelok;

		recursetrackpossible(individ* mother, const threadblock& tb, MarkerVal markerval, int marker, int upflag, int upflag2, int upshift, int genwidth, int f2n, int firstpar, int numrealpar, const trackpossibleparams& extparams) :
		genwidth(genwidth), extparams(extparams), markerval(markerval), marker(marker), tb(tb), firstpar(firstpar), mother(mother)

		{
			upflagr = upflagit(upflag, firstpar, genwidth);
			upflag2r = upflagit(upflag2, firstpar, genwidth >> (NUMGEN - NUMFLAG2GEN));
			upshiftr = upflagit(upshift, firstpar, genwidth >> (NUMGEN - NUMSHIFTGEN));

			prelok = true;
		}

		operator double()
		{
			if (!prelok)
			{
				return 0;
			}

			double baseval =
				mother->pars[firstpar]->trackpossible<update, zeropropagate>(tb, markerval, marker,
				upflagr,
				upflag2r,
				upshiftr, extparams, genwidth >> 1);


			return baseval;
		}
	};

	// zeropropagate also implies that there is a gstr value, we propagate zeros to find any possible source strain
	// The main logic of tracking a specific inheritance pathway, computing haplotype weights, and overall feasibility.
	// update: Should haplotype weights be updated?
	// zeropropagate: Should zero marker values be kept, or changed to the actual values they are matched against.
	// threadblock: Reference to all thread-local data.
	// inmarkerval: The marker val we are matching against.
	// marker: The marker number.
	// flag: The actual genotype flag. Note that this is zero-extended to be one bit more than the actual enumeration of flags
	// (lowest bit always 0).
	// flag99: What's mostly called flag2. A specific shift state of marker numbers, or -1 to indicate that all shifts are allowed.
	// localshift: Based on shiftflagmode, the overall mapping *for the complete sequence* of strand numbers to actual parentage for all
	// individuals in the analysis.
	// extparams: external parameters.
	// genwidth: The width of the generation flags.
	template<bool update, bool zeropropagate> double trackpossible(const threadblock& tb, MarkerVal inmarkerval, const unsigned int marker,
		const unsigned int flag, const int flag99, int localshift = 0, const trackpossibleparams& extparams = tpdefault,
		const int genwidth = 1 << (NUMGEN - 1)) /*const*/
	{
		if (this == NULL) return 1;

		int upflag2 = -1;
		const int upflag = flag >> 1;
		const int upshift = localshift >> 1;
		int f2s = 0;
		const MarkerVal* themarker = &markerdata[marker].first;
		float themarkersure[2] = {0, 0};
		bool allthesame = themarker[0] == themarker[1];
		int f2end = 2;
		if (flag99 != -1 && genwidth >> (NUMGEN - NUMFLAG2GEN) > 0)
		{
			upflag2 = flag99 >> 1;
			f2s = flag99;
			f2end = flag99 + 1;			
		}

		int firstpar = flag & 1;
		double ok = 0;


		// flag2 determines which value in the tuple to check. flag and localshift determine the haplotype weight value
		// assigned to that test, and which parent to match against that value.
		for (int flag2 = f2s; flag2 < f2end && (HAPLOTYPING || !ok); flag2++)
		{

			MarkerVal markerval = inmarkerval;
			double baseval;

			int f2n = (flag2 & 1);
			int realf2n = f2n;

			// If this marker value is not compatible, there is no point in trying.
			if (markermiss<zeropropagate>(markerval, themarker[f2n]))
			{
				baseval = themarkersure[f2n];
			}
			else
			{
				baseval = 1.0 - themarkersure[f2n];
			}

			if (!baseval) continue;

			// Normalize, in some sense.
			f2n ^= ((firstpar ^ localshift) & 1);

			if (zeropropagate || !genwidth)
			{
				baseval *= 0.5;
			}
			else if (allthesame && themarkersure[0] + themarkersure[1] == 0)
			{
				baseval *= ((f2n) ? 1.0 : 0.0);
			}
			else
			{
					// No haplotype weights, all interpretations allowed.
					baseval *= 0.5;
			}



			if (!baseval)
			{
				continue;
			}

			// There should be some other flag for the actual search depth
			if (genwidth == HAPLOTYPING)
			{
				if (zeropropagate && extparams.gstr)
				{
					*(extparams.gstr) *= 2;
					*(extparams.gstr) += strain - 1;
				}
			}
			else
			{
				// Track to the parent generation, creating evaluation objects first
				// These do a lookup in a special hash for combinations known to be 0, to avoid unnecessary calls
				// Both are checked before either branch of the pedigree is actually traced.
				recursetrackpossible<update, zeropropagate> subtrack1 =
					recursetrackpossible<update, zeropropagate>(this, tb, markerval, marker,
					upflag,
					upflag2,
					upshift,
					genwidth,
					f2n,
					firstpar,
					0,
					extparams);

				if (subtrack1.prelok && (!zeropropagate || (genwidth == 1 << (NUMGEN - 1))) )
					baseval *= recursetrackpossible<update, zeropropagate>(this, tb, themarker[!realf2n], marker,
					upflag,
					upflag2,
					upshift,
					genwidth,
					f2n,
					!firstpar,
					1,
					extparams);
				if (!baseval) continue;

				baseval *= subtrack1;
			}

			if (!baseval) continue;

			ok += baseval;

		}

		return ok;
	}


	// calltrackpossible is a slight wrapper that hides at least some of the interanl parameters needed for the recursion from outside callers
	template<bool update, bool zeropropagate> double calltrackpossible(const threadblock& tb, const MarkerVal* markervals, const unsigned int marker,
		const int genotype, const unsigned int offset, const int flag2, const double updateval = 0.0)
	{		
		return trackpossible<update, zeropropagate>(tb, UnknownMarkerVal,
			marker, genotype * 2, flag2, *(tb.shiftflagmode), trackpossibleparams(updateval, 0));
	}

	// "Fix" parents, i.e. infer correct marker values for any zero values existing.
	// Depending on some specifics, this can either assume that all existing marker values will be reflected in some offspring individuals,
	// or indeed only infer specifically those things that can be said from certain, i.e. if a zero value exists, and the other branch of
	// the pedigree doesn't allow a marker value found in the offspring, then it can be assumed to have been transmitted from the individual
	// with the zero value. 
	void fixparents(unsigned int marker, bool latephase)
	{
		MarkerValPair& themarker = markerdata[marker];
		individ* parp[3] = {pars[0], pars[1], this};
		int az = 0;

		for (int i = 0; i < 3; i++)
		{
			if (parp[i])
			{
				az += parp[i]->markerdata[marker].first == UnknownMarkerVal;
				az += parp[i]->markerdata[marker].second == UnknownMarkerVal;
				while (parp[i]->markervals.size() <= marker)
				{
					parp[i]->markervals.resize(markerdata.size());
				}
			}
		}


		double okvals[2] = {0};
		// We only accept an interpretation when it is by exclusion the only possible one. As soon as one intepretation has gained acceptance,
		// no need to retest it.
		for (shiftflagmode = 0; shiftflagmode < 1; shiftflagmode+=1)
		{
			threadblock tb;
			for (int i = 0; i < NUMTYPES; i++)
			{
				for (int flag2 = 0; flag2 < NUMPATHS; flag2++)
				{
					if (okvals[flag2 & 1]) continue;

					double ok = calltrackpossible<false, false>(tb, 0, marker, i, 0, flag2);
					if (!ok) continue;
					okvals[flag2 & 1] += ok;
					if (okvals[0] && okvals[1]) break;
				}
				if (okvals[0] && okvals[1]) break;
			}
		}

		if (!okvals[0] && !okvals[1])
		{
			printf("Clearing %d:%d\n", this->n, marker);
			markerdata[marker] = make_pair(UnknownMarkerVal, UnknownMarkerVal);
		}

		if ((((bool) okvals[0]) ^ ((bool) okvals[1])) || latephase)
		{
			for (int flag2 = 0; flag2 < 2; flag2++)
			{
				if (!okvals[flag2]) continue;

				for (int k = 0; k < 2; k++)
				{
					if (pars[k])
					{
						int u = ((k ^ flag2 /*^ *tb.shiftflagmode*/) & 1);
						if (latephase || (&themarker.first)[u] != UnknownMarkerVal)
						{
#pragma omp critical (parmarkerval)
							pars[k]->markervals[marker].insert((&themarker.first)[u]);
						}
					}
				}
			}
		}		
	}

	// Verifies that an update should take place and performs it. Just a wrapper to calltrackpossible.
	void updatehaplo(const threadblock& tb, const unsigned int marker, const unsigned int i, const int flag2, const double updateval)
	{
		MarkerValPair& themarker = markerdata[marker];		

		double ok = calltrackpossible<false, false>(tb, &themarker.first, marker, i, 0, flag2, updateval);

		if (ok)
		{
			calltrackpossible<true, false>(tb, &themarker.first, marker, i, 0, flag2, updateval);
		}
		else
		{

		}
	}

	// Adjust the probability, i.e. filter all probability values based on the haplotype weights and overall admissibility for the different
	// states.
	void adjustprobs(const threadblock& tb, PerStateArray<double>::T& probs, const unsigned int marker, double& factor, const bool oldruleout, int flag99)
	{
		double sum = 0;
		PerStateArray<double>::T probs2;

		const MarkerValPair& themarker = markerdata[marker];
		const bool ruleout = true; // TODO

		for (int q = 0; q <= (int) !ruleout; q++)
		{
			int f2start = 0;
			int f2end = NUMPATHS;

			// Negative values other than -1 just cause trouble

			// Can we really trust the logic to expand correctly even for zero values?
			//if (flag99 >= 0 || !somezero[marker])
			{
				f2start = flag99;
				f2end = flag99 + 1;
			}

			for (unsigned int i = 0; i < NUMTYPES; i++) // genotype
			{
				probs2[i] = probs[i];

				// We will multiply this already small number with an even smaller number... let's assume it's zero and be done with it.
				if (probs[i] < 1e-60)
				{
					probs[i] = 0;
					continue;
				}

				double realok = 0;

				for (int flag2 = f2start; flag2 < f2end; flag2++)
				{
					realok += calltrackpossible<false, false>(tb, &themarker.first, marker, i, 0, flag2);
				}					

					probs[i] *= (bool) realok;
				sum += probs[i];
			}

			if (sum == 0 && !ruleout)
			{
				for (int i = 0; i < NUMTYPES; i++)
				{
					probs[i] = probs2[i];
				}

				// The code sees: This doesn't make sense, ignore this marker!
				// NOTE: If the parent-grandparent genotypes are completely incompatible, the problem
				// is NOT eliminated.

				// Current caching scheme does not execute full elimination. Efficiency gains possible by not allowing
				// for genotyping error in every single marker. A small epsilon should otherwise be introduced whenever
				// any probs[i] is zero.
				markerdata[marker] = make_pair(UnknownMarkerVal, UnknownMarkerVal);
				fprintf(stderr, "Error in %x, marker %d, impossible path\n", this, marker);
				sum = 1;
			}
			else
				break;
		}

		// Normalize, update auxiliary exponent
		for (int i = 0; i < NUMTYPES; i++)
		{
			probs[i] /= sum;
		}
		factor += log(sum);
	}

	// Append a "multi-step" transition. If the cached values (essentially a N * N transition matrix for the steps from startmark to
	// endmark) are missing, calculate them first.
	double fillortake(const threadblock& tb, const int index, const unsigned int startmark, const unsigned int endmark, PerStateArray<double>::T& probs)
	{
		if ((tb.done[*(tb.shiftflagmode)])[index] != (*tb.generation))
		{

			for (unsigned int i = 0; i < NUMTYPES; i++)
			{
				PerStateArray<double>::T probs2 = {{0}};
				probs2[i] = 1;

				// Note ruleout here, could be set to false if we "prime" with an evenly distributed probs first
				(tb.factors[*tb.shiftflagmode])[index][i] = quickanalyze<false, noneturner>(tb, none, startmark,
					endmark,
					NONESTOP,
					-1,
					true,
					probs2);

				double sum = 0;

#pragma ivdep:back
				for (int j = 0; j < NUMTYPES; j++)
				{
					if (!_finite(probs2[j])) probs2[j] = 0.0;
					sum += probs2[j];
				}

				sum *= NUMTYPES;

				if (sum == 0) sum = 1;

				(tb.factors[*tb.shiftflagmode])[index][i] += log(sum);
				float* probdest = &(tb.cacheprobs[*tb.shiftflagmode])[index][i][0];

#pragma ivdep
				for (int j = 0; j < NUMTYPES; j++)
				{
					probdest[j] = probs2[j] / sum;
				}
			}
			(tb.done[*tb.shiftflagmode])[index] = (*tb.generation);
		}

		float factor = MINFACTOR;
		for (int i = 0; i < NUMTYPES; i++)
		{
			factor = max(factor, (tb.factors[*tb.shiftflagmode])[index][i]);
		}

		PerStateArray<double>::T probs2 = {{0}};

		for (int i = 0; i < NUMTYPES; i++)
		{
			float step = (tb.factors[*tb.shiftflagmode])[index][i] - factor;
			if (probs[i] == 0.0 || step < -100.0f) continue;
			double basef = exp((double) step) * probs[i];
			//			if (basef == 0.0 || !_finite(basef)) continue;
			const float* probsource = &(tb.cacheprobs[*tb.shiftflagmode])[index][i][0];
#pragma ivdep
			for (int j = 0; j < NUMTYPES; j++)
			{
				const double term = basef * probsource[j];
				probs2[j] += term;
			}
		}

		double sum = 0;

		for (int i = 0; i < NUMTYPES; i++)
		{
			sum += probs2[i];
		}

		if (sum == 0)
		{
			return MINFACTOR;
		}

		factor += log(sum);
		sum = 1 / sum;

		for (int i = 0; i < NUMTYPES; i++)
		{
			probs[i] = probs2[i] * sum;
		}

		return factor;
	}

	// Analyze for a specific range, including a possible fixed specific state at some position (determined by stopdata)
	template<bool inclusive, class T, class G> double quickanalyze(const threadblock& tb, const T& turner, unsigned int startmark,
		const unsigned int endmark, const G &stopdata, const int flag2, bool ruleout, PerStateArray<double>::T& probs,
		float minfactor = MINFACTOR)
	{
		unsigned int stepsize;
		double factor = 0;
		bool allowfull = inclusive;
		bool frommem = false;


		// Loop until we have reached the end, with blocks of varying sizes.
		while (startmark < endmark)
		{
			for (stepsize = 1; stepsize < (endmark - startmark + allowfull) &&
				stopdata.okstep(startmark, startmark + stepsize) &&
				!(startmark & (stepsize - 1)); stepsize *= 2);

			// A single step, either due to the fixated genotypee being within this range, or simply because we've gone all the way down
			// the tree.
			if (stepsize <= 2)
			{
				stepsize = 1;

				bool willquickend = false;

					factor += realanalyze<0, T>(tb, turner, startmark, startmark + stepsize, stopdata, flag2, ruleout, &probs);

				// This will work.
				if (!_finite(factor) || factor <= minfactor)
				{
					return MINFACTOR;
				}
			}
			else
			{
				// Use a stored multi-step transition.
				stepsize /= 2;

				int base = 0;
				for (int q = stepsize / 2; q > 1; q /= 2)
				{
					base += markerposes.size() / q;
				}

				base += startmark / stepsize;

				factor += fillortake(tb, base, startmark, startmark + stepsize, probs);
			}
			startmark += stepsize;
			allowfull |= true;

			if (!_finite(factor) || factor <= minfactor)
			{
				return MINFACTOR;
			}
		}		

		return factor;
	}

	// A wrapper to quickanalyze, preparing the start probability vector.
	template<class T, class G> double doanalyze(const threadblock& tb, const T& turner, const int startmark, const int endmark, const G& stopdata,
		const int flag2, bool ruleout = false, PerStateArray<double>::T* realprobs = 0, float minfactor = MINFACTOR)
	{
		PerStateArray<double>::T fakeprobs;
		PerStateArray<double>::T& probs = realprobs ? *realprobs : fakeprobs;

		if (realprobs)
		{
			//probs = realprobs;
		}
		else
		{
			//probs = fakeprobs;
			for (int i = 0; i < NUMTYPES; i++)
			{
				fakeprobs[i] = EVENGEN;
			}
		}

		double factor = quickanalyze<true, T>(tb, turner, startmark, endmark, stopdata, flag2, ruleout, probs, minfactor);
		bool small = !_finite(factor) || minfactor >= factor;

		if (!small) adjustprobs(tb, probs, endmark, factor, ruleout, -1); // TODO, the very last marker can be somewhat distorted

		return factor;
	}	

	// This is the actual analyzing code. It works with no caches, and can function independently, but is generally only used to patch in those
	// elements not found in caches by quickanalyze and fillortake.
	//
	// Both transition and emission (through adjustprobs) probabilities are handled here.
	//
	// first bit in updateend signals whether the interval is end-inclusive at endmark
	// the second bit in updateend will be 0 if the interval is end-inclusive at startmark, and 1 IF NOT
	// the third bit will cause the code to quit early, after processing the genotype and turner condition
	template<int updateend, class T, class G> double realanalyze(const threadblock& tb, const T& turner, const int startmark, const int endmark, const G& stopdata,
		const int flag2, const bool ruleout = false, PerStateArray<double>::T* realprobs = 0)
	{
		PerStateArray<double>::T fakeprobs;
		PerStateArray<double>::T& probs = realprobs ? *realprobs : fakeprobs;

		if (realprobs)
		{
			//probs = realprobs;
		}
		else
		{
			//probs = fakeprobs;
			for (int i = 0; i < NUMTYPES; i++)
			{
				fakeprobs[i] = EVENGEN;
			}
		}

		// The *logged* normalization factor
		double factor = 0;

		// Walk over all markers.
		for (int j = startmark + 1; j <= endmark; j++)
		{
			double startpos = markerposes[j - 1];
			double endpos = markerposes[j];
			int genotype = -1;

			bool tofind = stopdata.fixtofind(genotype, startpos, endpos, j);

			int f2use = -1;

			if (tofind)
			{
				f2use = flag2;
			}

			if (genotype != -2)
			{
				// If we are at the very first position, and the specific flag was set, include the emission probabilities for the previous
				// marker. Used to maximize the caching.
				if (!((updateend & 2) && (j == startmark + 1))) adjustprobs(tb, probs, j - 1, factor, ruleout, f2use);
			}
			else
			{
				// We could do some stuff here to avoid excessive adjustprobs calls
				// a -2 genotype does not only mean that all genotypes are allowed, but indeed that the marker data at this marker
				// is ignored!
			}

			// For a specific intra-marker region, we have two cases: the case of a fixated position between the two markers, and the simple case
			// of no fixated position.
			for (int iter = 0; iter <= (int) tofind; iter++)
			{
				// If iter is 1, we have currently handled the transition all the way to the fixated position. Now filter to keep only
				// a single state value positive.
				if (iter)
				{
					turner(probs);
					if (genotype >= 0)
					{
						for (int i = 0; i < NUMTYPES; i++)
						{
							probs[i] = probs[i] * (i == genotype);
						}
					}

					// Were we asked to stop at this very state, in the middle of things?
					if (updateend & 4) return factor;
				}


				double dist = (endpos - startpos);

				// Compute transitions, assuming there is any distance to transfer over.
				if (dist > 0)
				{
					PerStateArray<double>::T probs2 = {{0}};
					double recprob[2];

#pragma ivdep
					// Compute recombination probabilities for this specific distance, for the two sexes.
					// (as the sex-dependent marker distance might not be a simple transformation, the actrec
					// data comes into play).
					for (int k = 0; k < 2; k++)
					{
						// There is no recombination going on in the heterogametic 
						if (sexmarker[j] && k != homogamsex)
						{
							recprob[k] = 0;
						}
						else
						{						
							recprob[k] = 0.5 * (1.0 - exp(getactrec(stopdata, startpos, endpos, k, j) * (dist)));
							//					if (iter == tofind) recprob[k] = max(recprob[k], 1e-5);
						}
					}

					// Precompute values for presence (/lack of) a crossover event for either sex.
					double other[2][2];
#pragma ivdep
					for (int m = 0; m < 2; m++)
					{
						for (int k = 0; k < 2; k++)
						{
							double prob = recprob[k];
							if (m) prob = 1.0 - prob;

							other[m][k] = prob;
						}
					}

					PerStateArray<double>::T recombprec;
#pragma ivdep
					for (int index = 0; index < NUMTYPES; index++)
					{
						recombprec[index] = 1;
					}


					// Compute probabilities for arbitrary xor values of current and future state
					for (int t = 0; t < TYPEBITS; t++)
					{
						int sex = TYPESEXES[t];

#pragma ivdep
						for (int index = 0; index < NUMTYPES; index++)
						{
							int val = !((index >> t) & 1);
							recombprec[index] *= other[val][sex];
						}

					}

					// Use those xor values
					// For the 4-state model, this is an inefficient way to go about it, but it is quite a bit more efficient for
					// the 64-state model (or beyond).
					for (int from = 0; from < NUMTYPES; from++)
					{
						if (probs[from] < MINFACTOR) continue;
						for (int to = 0; to < NUMTYPES; to++)
						{
							probs2[to] += probs[from] * recombprec[from ^ to];
						}
					}

					for (int c = 0; c < NUMTYPES; c++)
					{
						probs[c] = probs2[c];
					}
				}
				//				else
				{
					//				if (iter == tofind)
					{
						for (int c = 0; c < NUMTYPES; c++)
						{
							if (probs[c] < 1e-20) probs[c] = 1e-20;
						}
					}
				}

				startpos = endpos;
				endpos = markerposes[j];
			}			
		}

		if (updateend & 1)
		{
			adjustprobs(tb, probs, endmark, factor, ruleout, -1); // TODO
		}

		return factor;
	}
};

// Oh, how we waste memory, in a pseudo-O(1) manner
individ* individer[1000000];
// dous contains those individuals that really should be analyzed
vector<individ*> dous;

// retrieve the individual with a specific number
individ* const getind(int n, bool real = false)
{
	if (n <= 0) return 0;	

	if (!real && !individer[n]) n = 0;

	if (!individer[n])
	{		
		individer[n] = new individ();
		individer[n]->n = n;
	}

	return individer[n];
}

miniindividual* const getmind(int n, bool real)
{
	return getind(n, real);
}


void initind(miniindividual* ind1)
{
	individ* ind = (individ*) ind1;

	ind->markerdata.resize(markerposes.size());
}


// Some operations performed when marker data has been read, independent of format.
void postmarkerdata()
{
	int any, anyrem;
	bool latephase = false;


	markerweight.resize(markerposes.size());
	// If inference is active, add "new" marker data until all data has been found.
	if (CORRECTIONINFERENCE) do
	{
		for (int i = 1; i < 1000000; i++)
		{
			individ* ind = getind(i);
			if (ind->markerdata.size())
			{
				generation++;

				for (int g = 0; g < ind->markerdata.size(); g++)
				{					
					ind->fixparents(g, latephase);					
				}
			}
		}

		any = 0;
		anyrem = 0;
		for (int i = 1; i < 1000000; i++)
		{
			individ* ind = getind(i);

			for (int g = 0; g < (int) ind->markervals.size(); g++)
			{
				if (g == 1) continue;

				int startsize = ind->markervals[g].size();
				const int known =
					(ind->markerdata[g].first != UnknownMarkerVal) +
					(ind->markerdata[g].second != UnknownMarkerVal);

				const int oldany = any;

				if (latephase && known == 2 && !startsize && ind->gen < 2)
				{
					ind->markerdata[g].second = UnknownMarkerVal;
					ind->markerdata[g].first = UnknownMarkerVal;
					any++;
					anyrem++;
				}
				else
					if (known == 2) continue;

				ind->markervals[g].erase(UnknownMarkerVal);
				startsize = ind->markervals[g].size();

				if (ind->markerdata[g].first  != UnknownMarkerVal) ind->markervals[g].insert(ind->markerdata[g].first);
				if (ind->markerdata[g].second != UnknownMarkerVal) ind->markervals[g].insert(ind->markerdata[g].second);

				if (ind->markervals[g].size() >= 3)
				{
					fprintf(stderr, "Error, too many matches: %d\t%d\n", i, g);
				}
				if (!latephase && ind->markervals[g].size() == 2)
				{
					ind->markerdata[g] = make_pair(*ind->markervals[g].begin(), *(++ind->markervals[g].begin()));
					any++;
				}
				if (latephase && ind->markervals[g].size() == 1 && startsize == 1 && known == 1 && !ind->pars[0] && !ind->pars[1])
				{
					if (ind->markerdata[g].first == UnknownMarkerVal || ind->markerdata[g].second == UnknownMarkerVal) any++;
					ind->markerdata[g] = make_pair(*ind->markervals[g].begin(), *ind->markervals[g].begin());					
				} // DANGEROUS ASSUMPTIONS
				else if (!latephase && ind->markervals[g].size() == 1 && known == 0)
				{
					any++;
					ind->markerdata[g] = make_pair(*ind->markervals[g].begin(), UnknownMarkerVal);
				}

				
				if (any != oldany && ind->markerdata[g].first.value()) printf("Correction at %d, marker %d (%d;%d)\n", i, g,
					ind->markerdata[g].first.value(), ind->markerdata[g].second.value());

			}

			for (int g = 0; g < (int) ind->markervals.size(); g++)
			{
				if (ind->markerdata[g].first == sexmarkerval) {
					ind->markerdata[g] = make_pair(ind->markerdata[g].second, ind->markerdata[g].first);
				}
				ind->markervals[g].clear();
			}
		}
		fprintf(stderr, "Number of corrected genotypes: %d\n", any);
		if (latephase)
		{
			latephase = false;
		}
		else
		{
			if (!any && !latephase)
			{
				any++;
				latephase = true;
			}
		}
	}
	while (any > anyrem);

	for (int i = 1; i < 1000000; i++)
	{
		individ* ind = getind(i);
		ind->markervals.clear();

	}
}


void add_dous(miniindividual* ind)
{
	dous.push_back((individ*) ind);

}

void clear_dous()
{
	for (int i = 0; i < 1000000; i++)
	{
		individ* ind = individer[i];
		if (ind)
		{
			delete ind;
			individer[i] = 0;
		}
	}
	dous.clear();
}

// The actual walking over all chromosomes for all individuals in "dous"
// If "full" is set to false, we assume that haplotype inference should be done, over marker positions.
// A full scan is thus not the iteration that takes the most time, but the scan that goes over the full genome grid, not only
// marker positions.
template<bool full> void doit(resultreceiver& receiver
#ifdef F2MPI
	, mpi::communicator& world
#endif
	)
{
	const bool doprint = full;
#ifdef F2MPI
	broadcast(world, actrec, 2, 0);
#endif

	int count = 0;
	vector<vector<array<float, 2> > > realgeno;

	realgeno.resize(dous.size());

	for (unsigned int j = 0; j < dous.size(); j++)
	{
		if (dous[j]->markerdata.size())
		{
			count++;
		}
		else
		{
			fprintf(stderr, "ZERO RESULTS %d\n", j);
			dous.erase(dous.begin() + j);
		}
	}

	static int iter = 0;
	iter++;
	/*if (iter == 35)
	{
	individ* inds[2] = {getind(2348), getind(2352)};
	for (int j = 0; j < 2; j++)
	{
	individ* ind = inds[j];
	for (int i = 4; i < ind->haploweight.size(); i++)
	{
	ind->haploweight[i] = 0.5;
	}
	}
	}*/

	for (int i = 0; i < 1000000; i++)
	{
		individ* ind = getind(i);
		if (!ind) continue;
		ind->children = 0;
		ind->kids.clear();
		//		fflush(out);

	}

	for (int j = 0; j < (int) dous.size(); j++)
	{

		for (int i = 0; i < 2; i++)
		{
			individ* pnow = dous[j]->pars[i];

			if (pnow)
			{
				pnow->children++;
				pnow->kids.push_back(dous[j]);
			}
		}
	}


	for (unsigned int i = 0; i < chromstarts.size() - 1; i++)
	{
		for (int j = 0; j < (int) dous.size(); j++)
		{			
#ifdef F2MPI
			if (j % world.size() != world.rank()) continue;			
#endif
			//			realgeno[j].resize(markerposes[chromstarts[1] - 1] + 1);

			generation++;
			threadblock tborig;
			threadblock tb = tborig;

			for (int t = 0; t < NUMSHIFTS; t++)
			{
				factors[t].resize(markerposes.size());
				cacheprobs[t].resize(markerposes.size());
				done[t].resize(markerposes.size());
			}

			if (dous[j]->markerdata.size())
			{				
				int qstart = -1000 - chromstarts[i];
				int qend = -1000 - chromstarts[i + 1];
				int qd = -1;
				int f2s = 0;
				int f2end = NUMPATHS;

					f2s = -1;
					f2end = 0;

				int shifts = 0;
				int shiftend = NUMSHIFTS;

				reltree.clear();
				reltree.push_back(dous[j]);
				int flag2ignore = 0;


				sort(reltree.begin(), reltree.end());
				reltree.resize(unique(reltree.begin(), reltree.end()) - reltree.begin());

				bool skipsome = false;
				for (int u = 0; u < reltree.size() && !skipsome; u++)
				{
					for (int p = 0; p < 2; p++)
					{
						if ((&(reltree[u]->markerdata[chromstarts[i]].first))[p] == sexmarkerval)
							skipsome = true;
					}
				}

					reltree.resize(0);
					flag2ignore = 0;

				if (full)
				{
					qstart = (int) markerposes[chromstarts[i]];
					qend = (int) markerposes[chromstarts[i + 1] - 1] + 1;
					qd = 1;
					f2s = -1;
					f2end = 0;
					/*shifts = 0;
					shiftend = 1;*/
				}			

				if (dous[j]->gen < 2) shiftend = min(2, shiftend);

				double factor = -1e15;
				double factors[NUMSHIFTS];
				for (shiftflagmode = shifts; shiftflagmode < shiftend; shiftflagmode++)
				{
					factors[shiftflagmode] = dous[j]->doanalyze<noneturner>(tb, none, chromstarts[i], chromstarts[i + 1] - 1, NONESTOP, -1, false, 0, -100 + factor);
					factor = max(factor, factors[shiftflagmode]);
				}

				// Normalize!
				double realfactor = 0;
				for (int s = shifts; s < shiftend; s++)
				{
					realfactor += exp(factors[s] - factor);
				}


				factor += log(realfactor);
				if (_isnan(factor)) continue;

				// States are mapped onto values describing the line/strain origin, in the sense of 00, 01, 10 or 11
				PerStateArray<int>::T maptogeno;
				shiftflagmode = 0;
				for (int g = 0; g < NUMTYPES; g++)
				{
					int sum = 0;
					dous[j]->trackpossible<false, true>(tb, UnknownMarkerVal
						, 0, g * 2, 0, 0, trackpossibleparams(0, &sum));

					// switch due to earlier inner switching
					if (sum == 2 || sum == 1) sum = 3 - sum;

					maptogeno[g] = sum;
				}

				// Walk over all chromosome positions, whether it be markers (negative q values <= -1000) or grid positions


				for (int q = qstart; q != qend; q+=qd)
				{
					double probs[4] = {0};
					double mwvals[NUMTYPES][NUMTYPES] = {0};
					double mwfvals[NUMTYPES] = {0};

					double mwval[4] = {0};

					for (int g = 0; g < NUMTYPES; g++)
					{						
						for (shiftflagmode = shifts; shiftflagmode < shiftend; shiftflagmode++)
						{
							if (factor - factors[shiftflagmode] > 10) continue;

							for (int flag2 = f2s; flag2 < f2end; flag2++)
							{
								if (flag2 & (flag2ignore)) continue;

								int firstpar = 0;
								double val;


								val = dous[j]->doanalyze<noneturner>(tb, none, chromstarts[i], chromstarts[i + 1] - 1, classicstop(q, g),
									flag2, true, 0, -100.0 + factor) - factor;

								if (_finite(val) && val > -20.0)
								{
									int f2n = ((flag2/* ^ shiftflagmode*/) & 1);
									val = exp(val);


									int mapval = maptogeno[g];
									// int g3 = (g & 1) + ((bool) (g & 8)) * 2;
									receiver.pushresult(dous[j], i, q, mapval, g, flag2, shiftflagmode, val);
								}
continueloop:;
							}
						}
					}



				}
			}
		}

	}



	generation++;
}

// Handle the fuss of trailing \r, \n characters when combining scanf and gets and stuff.
void clean(char* tlf)
{
	int i = strlen(tlf);
	while (i && tlf[i - 1] < 32)
	{
		tlf[i - 1] = 0;
		i--;
	}
}



template void doit<true>(resultreceiver&);
