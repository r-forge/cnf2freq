#include <vector>
#include <boost/units/quantity.hpp>
#include <boost/units/conversion.hpp>
#include <boost/units/make_system.hpp>
#include <boost/units/base_unit.hpp>
#include <boost/units/base_dimension.hpp>
#include <boost/units/static_constant.hpp>

#include <boost/integer_traits.hpp>
// #define MCIBD

using namespace boost;
//using namespace boost::mpi;
using namespace boost::units;

// Declarations used for type safety based on Boost Units
class FactorUnit {};

extern bool early;
struct MarkerBaseDimension : base_dimension<MarkerBaseDimension, 931> {};
typedef MarkerBaseDimension::dimension_type MarkerDimension;
struct MarkerValBaseUnit : public base_unit<MarkerValBaseUnit, MarkerDimension, 932> {};
typedef make_system<MarkerValBaseUnit>::type UnitSystem;

typedef unit<MarkerDimension, UnitSystem> MarkerValUnit;

BOOST_UNITS_STATIC_CONSTANT(MarkerValue,MarkerValUnit);  


typedef quantity<MarkerValUnit, int> MarkerVal;
typedef quantity<FactorUnit, float> Factor;

typedef std::pair<MarkerVal, MarkerVal> MarkerValPair;

const MarkerVal UnknownMarkerVal = (MarkerVal) 0;
const MarkerVal sexmarkerval = (boost::integer_traits<MarkerVal::value_type>::const_max - 1) * MarkerValue;


struct miniindividual
{
	// The individual #.
	int n;
	// Generation number. Convenient, while not strictly needed.
	int gen;	
	// Parents.
	miniindividual* minipars[2];
	
	// Sex.
	bool sex;
	// Line or strain of origin, should only exist in founders.
	int strain;
	// Marker data as a list of pairs. No specific ordering assumed.
	std::vector<MarkerValPair > markerdata;	
	// The haplotype weight, or skewness. Introducing an actual ordering of the value in markerdata.
	std::vector<float> haploweight;

	std::vector<std::pair<float, float> > markersure;	
};

void postmarkerdata();
miniindividual* const getmind(int n, bool real = false);

extern int homogamsex;
extern std::vector<double> markerposes;
extern std::vector<bool> sexmarker;
extern std::vector<double> actrec[2];
extern std::vector<unsigned int> chromstarts;

// Note that all calls to this method can be done in a reentrant manner
// Read: "can" means "implementations need to support this to work in a multithreaded version"
// Read: If MPI is used, different instances will actually be used and boost serialization will
// take place. Whoah!
struct resultreceiver
{
	virtual void pushresult(miniindividual* ind, int chromo, int q, int mapval, int g, int flag2, int shiftflagmode, double val);
};

void add_dous(miniindividual*);
void clear_dous();

template<bool full> void doit(resultreceiver& receiver
#ifdef F2MPI
	, mpi::communicator& world
#endif
	);

void initind(miniindividual* ind1);