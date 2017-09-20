/**** Utility class useful for working with BioC structures
 ****/

/**** Node_Converter
      Copies data elements from old object to new while descending the
      tree.

      Override methods where something interesting should happen.
****/

#include "BioC.hpp"

namespace BioC {

  class Node_Converter {
  public:

    virtual void convert( const Relation & relation,
                          Relation & new_relation );
    virtual void convert( const Annotation & annotation,
                          Annotation & new_annotation );
    virtual void convert( const Sentence & sentence,
                          Sentence & new_sentence );
    virtual void convert( const Passage & passage,
                          Passage & new_passage );
    virtual void convert( const Document & document,
                          Document & new_document );
    virtual void convert( const Collection & collection,
                          Collection & new_collection );
  };

}
