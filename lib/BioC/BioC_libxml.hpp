
/****  Class to read / write IRET Collection classes using libxml
 ****/

#include <string>

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <libxml/xmlreader.h>
#include <libxml/xmlwriter.h>

#include "BioC.hpp"

using std::map;
using std::string;

namespace BioC {

  static const string default_dtd = "BioC.dtd";

  class Connector_libxml {
  public:
    Connector_libxml(void) :
    dtd(default_dtd)
    {}

    string dtd;                 // set to default in constructor
    // defined by reading; set before writing, if needed

    // in memory

    void read(const string & docname, Collection & collection );
    void write(const string & docname, const Collection & collection );

    // serial

    void start_read(const string & docname, Collection & collection );
    bool read_next( Document & document );
    void start_write( const string & docname, const Collection & collection );
    void write_next( const Document & document );
    void end_write( void );

    // internal in memory methods

    void from_XML( xmlNodePtr cur, Collection & collection );
    void from_XML( xmlNodePtr cur, Document & document );
    void from_XML( xmlNodePtr cur, Passage & passage );
    void from_XML( xmlNodePtr cur, Sentence & sentence );
    void from_XML( xmlNodePtr cur, Annotation & annotation );
    void from_XML( xmlNodePtr cur, Relation & relation );

    // return true if node present
    bool getValue( xmlNodePtr cur, const char * element, string & value );
    bool getValue( xmlNodePtr cur, const char * element, int & value );
    bool getValue( xmlNodePtr cur, const char * element,
                   vector<string> & value );
    bool getValue( xmlNodePtr cur, const char * element,
                   map<string,string> & value );
    bool getValue( xmlNodePtr cur, const char * element,
                   vector<Location> & value );
    bool getValue( xmlNodePtr cur, const char * element,
                   vector<Node> & value );
    bool getAttribute( xmlNodePtr cur, const string & name, string & value );
    bool getAttribute( xmlNodePtr cur, const string & name, int & value );


    void to_XML( const Collection & collection );
    void to_XML( xmlNodePtr root, const Document & document );
    void to_XML( xmlNodePtr document_xml, const Passage & passage );
    void to_XML( xmlNodePtr passage_xml, const Sentence & sentence );
    void to_XML( xmlNodePtr sentence_xml, const Annotation & annotation );
    void to_XML( xmlNodePtr annotation_xml, const Location & location );
    void to_XML( xmlNodePtr sentence_xml, const Relation & relation );
    void to_XML( xmlNodePtr relation_xml, const Node & node );

    template< typename Type >
    xmlNodePtr setValue( xmlNodePtr parent, const string & name,
                         const Type & value ) {
      if ( value.size() <= 0 )
        return 0;
    
      for ( int i = 0; i < value.size(); ++i )
        to_XML( parent, value[i] );
      return 0;
    }
    void setAttribute( xmlNodePtr xml_node, const string & name,
                       const string & value );
    void setAttribute( xmlNodePtr xml_node, const string & name,
                       int value );

    void parse_error( int return_code );

    // internal serial

    int getDocument( Document & document );
    int getValue( const string & element, int & value );
    int getValue( const string & element, string & value );
//    string readAttribute( xmlNodePtr cur, const string & name );
    
    void write_next( const Passage & passage );
    void write_next( const Sentence & sentence );
    void write_next( const Annotation & annotation );
    void write_next( const Location & location );
    void write_next( const Relation & relation );
    void write_next( const Node & node );

    void endElement( void );
    void startElement( const string & name );
    void writeAttribute( const string & name, const string & value );
    void writeAttribute( const string & name, int value );
    void writeElement( const string & name, int value );
    void writeElement( const string & name, const string & value );
    void writeElement( const string & name, const vector<string> & value );
    void writeElement( const string & name, const map<string,string> & value );

    xmlDocPtr xml_doc;

    xmlTextReaderPtr reader;
    xmlTextWriterPtr writer;
  };

}
