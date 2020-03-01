// 2010 © Václav Šmilauer <eudoxos@arcig.cz>

#include <sudodem/pkg/common/MatchMaker.hpp>

SUDODEM_PLUGIN((MatchMaker));

Real MatchMaker::operator()(int id1, int id2, Real val1, Real val2) const {
	FOREACH(const Vector3r& m, matches){
		if(((int)m[0]==id1 && (int)m[1]==id2) || ((int)m[0]==id2 && (int)m[1]==id1)) return m[2];
	}
	// no match
	if(fbNeedsValues && (isnan(val1) || isnan(val2))) throw std::invalid_argument("MatchMaker: no match for ("+boost::lexical_cast<string>(id1)+","+boost::lexical_cast<string>(id2)+"), and values required for algo computation '"+algo+"' not specified.");
	return computeFallback(val1,val2);
}

void MatchMaker::postLoad(MatchMaker&){
	if(algo=="val")      { fbPtr=&MatchMaker::fbVal; fbNeedsValues=false; }
	else if(algo=="zero"){ fbPtr=&MatchMaker::fbZero;fbNeedsValues=false; }
	else if(algo=="avg") { fbPtr=&MatchMaker::fbAvg; fbNeedsValues=true;  }
	else if(algo=="min") { fbPtr=&MatchMaker::fbMin; fbNeedsValues=true;  }
	else if(algo=="max") { fbPtr=&MatchMaker::fbMax; fbNeedsValues=true;  }
	else if(algo=="harmAvg") { fbPtr=&MatchMaker::fbHarmAvg; fbNeedsValues=true; }
	else throw std::invalid_argument("MatchMaker:: algo '"+algo+"' not recognized (possible values: val, avg, min, max, harmAvg).");
}

Real MatchMaker::computeFallback(Real v1, Real v2) const { return (this->*MatchMaker::fbPtr)(v1,v2); }
