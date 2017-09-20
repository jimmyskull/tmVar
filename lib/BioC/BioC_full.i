%module BioC_full
%{
/* Includes the header in the wrapper code */
#define SWIG_FILE_WITH_INIT
#include "BioC.hpp"
#include "BioC_libxml.hpp"
#include "BioC_util.hpp"
%}
 
/* Parse the header file to generate wrappers */
%include "std_map.i"
%include "std_vector.i"
%include "std_string.i"

/*%extend std::basic_string<char> {
    const std::string substr(size_t pos) {
        std::string res;
        res = $self->substr(pos, 1);
        return res;
    } 
}*/

%include "BioC.hpp"
%include "BioC_libxml.hpp"
%include "BioC_util.hpp"

/*namesapce BioC {*/
%template(nodV) std::vector<BioC::Node>;
%template(locV) std::vector<BioC::Location>;
%template(stcV) std::vector<BioC::Sentence>;
%template(antV) std::vector<BioC::Annotation>;
%template(strV) std::vector<std::string>;
%template(psgV) std::vector<BioC::Passage>;
%template(rltV) std::vector<BioC::Relation>;
%template(docV) std::vector<BioC::Document>;
%template(strM) std::map<string,string>;

%apply const std::string& {std::string* };
/*namesapce BioC {
struct stdstring2perl
{
  std::string content;
}; 
}
/*}*/

