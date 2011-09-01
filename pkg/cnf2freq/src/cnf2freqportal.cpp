#define STRICT_R_HEADERS

// #include <boost/lexical_cast.hpp>
#include <Rcpp.h>
extern "C"
{
#include <R.h>
#include <Rdefines.h>
}
#include <cmath>

#include <map>
#include <vector>
#include <algorithm>

#include "cnf2freqmini.h"

const char* ROW_NAMES = "row.names";
const char* DIM_NAMES = "dimnames";

extern "C"
{
	SEXP cnf2doer(SEXP args);
}

const double ministep = 1E-5;

struct indmapper
{
	int now;
	typedef int eltype;
	typedef std::map<eltype, int> maptype;
	typedef std::map<int, eltype> maptypeback;

	maptype inds;
	maptypeback back;

	indmapper() : now(1)
	{			
	}

	int getrev(miniindividual* ind)
	{
		return back[ind->n];
	}

	miniindividual* getind(int whatter)
	{
		if (Rcpp::traits::is_na<INTSXP>(whatter))
		{
			return 0;
		}

		Rcpp::RObject r;

		maptype::iterator i = inds.find(whatter);
		int val;
		if (i == inds.end())
		{
			inds[whatter] = now;			
			val = now++;
			initind(getmind(val, true));
			back[getmind(val)->n] = whatter;
		}
		else
		{
			val = i->second;
		}
		
		return getmind(val);
	}

	/*miniindividual* getind(SEXP what)
	{					
		eltype whatter = Rcpp::as<eltype>(what);
		return getind(whatter);		
	}*/
};

struct genomapper
{
	int now;
	typedef int eltype;
	typedef std::map<eltype, int> maptype;

	maptype inds;

	genomapper() : now(1)
	{			
	}

	int getgene(int whatter)
	{					
		if (Rcpp::traits::is_na<INTSXP>(whatter))
		{
			return 0;
		}

		Rcpp::RObject r;

		maptype::iterator i = inds.find(whatter);
		int val;
		if (i == inds.end())
		{
			inds[whatter] = now;
			val = now++;
		}
		else
		{
			val = i->second;
		}
		
		return val;
	}
};

using namespace Rcpp;
struct cnf2doercontext
{	
	typedef std::vector<std::pair<std::pair<int, double>, int> > markerlist;
	indmapper indmap;
	genomapper genmap;

	GenericVector genodata;
	GenericVector phenodata;
	GenericVector data;

	std::vector<miniindividual*> our_dous;

	cnf2doercontext(GenericVector data) :
	data(data),
		genodata((SEXP) data["geno"]),
		phenodata((SEXP) data["pheno"])
	{
	}

	void createpedigree()
	{
		StringVector ids = phenodata.attr(ROW_NAMES);
		IntegerVector sex = phenodata["sex"];
		IntegerVector generation = phenodata["generation"];
		IntegerVector line = phenodata["line"];

		// Sex indicated here might not match actual sex, "sire" is parent 1 is parent 1
		// is parent 1.
		// (And as we are zero-indexed internally, that's actually parent 0.)
		IntegerVector idsire = phenodata["parent_1"];
		IntegerVector iddam = phenodata["parent_2"];

		for (int i = 0; i < ids.size(); i++)
		{
			int id;
			if (sscanf(ids[i], "%d", &id) != 1)
			{
			  //				Rf_error("Non-integer individual ID in geno matrix.");
			}

			miniindividual* ind = indmap.getind(id);
			miniindividual* sire;
			miniindividual* dam;


			sire = indmap.getind(idsire[i]);
			dam = indmap.getind(iddam[i]);

			ind->minipars[0] = sire;
			ind->minipars[1] = dam;
			ind->sex = sex[i] - 1;
			ind->strain = line[i];
			ind->gen = generation[i];
		}
	}

	void parsechromos(markerlist& markers, int heterogamsex, int sexchrom)
	{
		std::map<int, int> remap;
		for (int i = 0; i < markers.size(); i++)
		{
			remap[markers[i].second] = i;
		}

		// Walk over columns, because that's the way it's done
		std::set<miniindividual*> done;
		int endoffset = genodata.offset("chr");
		StringVector names = genodata.names();

		for (int i = 0; i < endoffset; i++)
		{
			int id;
			if (sscanf((char*) names[i], "%d", &id) != 1)
			{
			  //				Rf_error("Non-integer individual ID in geno matrix.");
			}

			miniindividual* ind = indmap.getind(id);
			int index = 0;
			if (done.find(ind) != done.end())
			{
				index++;
			}
			else
			{
				if (ind->gen == 3)
				{
					add_dous(ind);
					our_dous.push_back(ind);
				}
				done.insert(ind);
			}

			IntegerVector genotypes = genodata[i];
			for (int j = 0; j < markers.size(); j++)
			{
				(&(ind->markerdata[j].first))[index] = genmap.getgene(genotypes[markers[j].second]) * MarkerValue;

				if (index == 1 && ind->sex == heterogamsex && markers[j].first.first == sexchrom)
				{
					for (int k = 0; k < 2; k++)
					{
						MarkerVal& mv = (&(ind->markerdata[j].first))[(heterogamsex + k) & 1];
						if (mv == UnknownMarkerVal)
						{
							mv = sexmarkerval;
							break;
						}
					}					
				}
			}
		}
	}
};

void resultreceiver::pushresult(miniindividual* ind, int chromo, int q, int mapval, int g, int flag2, int shiftflagmode, double val)
{

}


// This implementation is not thread-safe in practice right now
struct lineresultreceiver : resultreceiver
{
	GenericVector& m;
	std::map<miniindividual*, int> linenos;
	const int size;

	lineresultreceiver(GenericVector& m, cnf2doercontext& ctxt, int size) : m(m), size(size)
	{
		for (int i = 0; i < ctxt.our_dous.size(); i++)
		{
			linenos[ctxt.our_dous[i]] = i * size;
		}
	}

	void pushresult(miniindividual* ind, int chromo, int q, int mapval, int g, int flag2, int shiftflagmode, double val)
	{
		NumericMatrix n = m[chromo];
#ifdef MCIBD
		n(linenos[ind] + g, q) += val;
#else
		n(linenos[ind] + mapval, q) += val;
#endif
	}

	void normalizecleanup()
	{
		for (int c = 0; c < m.size(); c++)
		{
			NumericMatrix n = m[c];

			for (int k = 0; k < n.ncol(); k++)
			{
				for (int i = 0; i < linenos.size(); i++)
				{
					double sum = 0;
					for (int j = 0; j < size; j++)
					{
						sum += n(i * size + j, k);
					}

					if (sum < 1e-10) sum = 1;
					sum = 1 / sum;
					
					for (int j = 0; j < size; j++)
					{
						n(i * size + j, k) *= sum;
					}
				}
			}
		}
	}
};

const char* F2suffixes[4] = {"11", "12", "21", "22"};
const int F2states = 4;
const int MCstates = 64;

SEXP createOutputVector(cnf2doercontext& ctxt, int length, const std::vector<std::string>& suffixes)
{
	NumericMatrix matr(ctxt.our_dous.size() * suffixes.size(), length);
	StringVector newRowNames(ctxt.our_dous.size() * suffixes.size());

	// The row names could be reused between chromos
	for (int i = 0; i < ctxt.our_dous.size(); i++)
	{
		for (int j = 0; j < suffixes.size(); j++)
		{
		  char tlf[255];
		  sprintf(tlf, "%d", ctxt.indmap.getrev(ctxt.our_dous[i]));
		  //std::string str = lexical_cast<std::string>(ctxt.indmap.getrev(ctxt.our_dous[i])) + std::string(".") + suffixes[j];
		  std::string str = std::string(tlf) + "." + suffixes[j];
		  newRowNames[i * suffixes.size() + j] = str;
		}
	}
	List dimNames(1);
	dimNames[0] = newRowNames;

	matr.attr(DIM_NAMES) = dimNames;

	return matr;
}

SEXP cnf2doer(SEXP args)
{
	using namespace Rcpp;

	clear_dous();
	markerposes.clear();
	sexmarker.clear();
	actrec[0].clear();
	actrec[1].clear();
	chromstarts.clear();

	GenericVector argv = args;
	GenericVector data = argv[1];
	GenericVector genodata = data["geno"];

	// Male and female monikers arbitrary, might no match actual data
	NumericVector malecm;
	NumericVector femalecm;
	NumericVector rulecm;
	IntegerVector chr  = genodata["chr"];
	int heterogam = as<int>(data["heterogam"]) - 1;	
	indmapper indmap;
	int sexchrom = as<int>(data["sex.chrom"]) - 1;

	double grid = 1.0;
	try
	{
	  grid = as<double>(argv["resolution"]);
	} catch (const index_out_of_bounds& ex){
		// Use default
	}
	
	malecm.set_sexp(genodata["sex_1_cM"]);
	femalecm.set_sexp(genodata["sex_2_cM"]);
	rulecm.set_sexp(genodata["ref_cM"]);

	cnf2doercontext::markerlist markers;
	// Repeated function calls. Consider global optimization/include-based introduction of all Rcpp components.
	// No micro-optimizations made to avoid cross-module boundary trivial function calls.
	for (int i = 0; i < malecm.size(); i++)
	{
		markers.push_back(std::make_pair(std::make_pair(chr[i] - 1, rulecm[i]), i));
	}

	// A very natural sorting
	std::sort(markers.begin(), markers.end());
	int oldChromo = -1;
	double last[3];

	// It makes more sense for the internal code to keep track of the homogametic sex.
	homogamsex = (heterogam == 1) ? 0 : 1;

	double invgrid = 1 / grid;

	for (int i = 0; i < markers.size(); i++)
	{
		if (oldChromo != markers[i].first.first)
		{
			oldChromo = markers[i].first.first;
			chromstarts.push_back(i);
			last[0] = 0;
			last[1] = 0;
			last[2] = -1;
		}

		int index = markers[i].second;
		double nowRule = rulecm[index] * invgrid;
		nowRule = std::max(last[2] + ministep, nowRule);

		markerposes.push_back(nowRule);
		sexmarker.push_back(markers[i].first.first == sexchrom);

		for (int j = 0; j < 2; j++)
		{
			if (i > 0)
			{
				actrec[0].push_back(std::min(-  (malecm[index] * invgrid - last[0]) / (nowRule - last[2]) * 0.02 * grid, -1e-5));
				actrec[1].push_back(std::min(-(femalecm[index] * invgrid - last[1]) / (nowRule - last[2]) * 0.02 * grid, -1e-5));
			}
			else
			{
				actrec[0].push_back(0);
				actrec[1].push_back(0);
			}
		}
		
		last[0] = malecm[index] * invgrid;
		last[1] = femalecm[index] * invgrid;
		last[2] = nowRule;
	}
	chromstarts.push_back(markers.size());

	cnf2doercontext ctxt(data);

	ctxt.createpedigree();
	ctxt.parsechromos(markers, heterogam, sexchrom);
	
	early = true;
	postmarkerdata();
	early = false;

	GenericVector outputvector(chromstarts.size() - 1);
	std::vector<std::string> suffixes;
#ifdef MCIBD
	suffixes.resize(MCstates);
	for (int i = 0; i < MCstates; i++)
	{
		  char tlf[255];
		  sprintf(tlf, "%d", i);
		  suffixes[i] = tlf;
	  //		suffixes[i] = lexical_cast<std::string>(i);
	}
#else
	suffixes.resize(F2states);
	for (int i = 0; i < F2states; i++)
	{
		suffixes[i] = F2suffixes[i];
	}
#endif
	for (int i = 0; i < chromstarts.size() - 1; i++)
	{
		outputvector[i] = createOutputVector(ctxt, (int) (markerposes[chromstarts[i + 1] - 1] + 1), suffixes);
	}

	lineresultreceiver x(outputvector, ctxt, suffixes.size());
	doit<true>(x);
	x.normalizecleanup();

	return outputvector;
}
