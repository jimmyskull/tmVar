#ifndef BIOC

/****  Classes for Collections of text Documents
 ****/

#include <iostream>
#include <map>
#include <string>
#include <vector>

using std::cout;
using std::map;
using std::ostream;
using std::string;
using std::vector;

namespace BioC {

  static const string bioc_dtd = "BioC.dtd";

  void write( ostream & out, const map<string,string> & infons );

  class Node {
  public:
    Node() {}
    Node( const string & in_refid, const string & in_role = "" ) :
      refid(in_refid), role(in_role)
    {}

    string refid;                 // id of Relation or Annotation
    string role;

    void write ( ostream & out = cout ) const;
  };

  ostream & operator<<( ostream & out, const Node & node );


  class Relation {
  public:
    string id;
    std::map<string,string> infons;
    std::vector<Node> nodes;

    int num_nodes(void) const { return nodes.size(); }
    const string & refid( int segment ) const { return nodes[segment].refid; }
    const string & role( int segment ) const { return nodes[segment].role; }
    void add_node( const string & refid, const string & role ) {
      add_node( Node(refid,role) );
    }
//    void add_node( Node & node ) { nodes.push_back(node); }
    void add_node( Node node ) { nodes.push_back(node); }

    void write ( ostream & out = cout ) const;
    void clear(void);
  };


  class Location {
  public:
    Location() {}
    Location( int in_offset, int in_length ) :
      offset(in_offset), length(in_length)
    {}

    int offset;
    int length;

    void write ( ostream & out = cout ) const;
  };

  ostream & operator<<( ostream & out, const Location & location );


  class Annotation {
  public:
  
    string id;
    std::map<string,string> infons;
    std::vector<Location> locations;
    string text;

    int num_segments(void) const;
    int offset( int segment = 0 ) const;
    int length( int segment = 0 ) const;
    void add_location( int offset, int length ) {
      add_location( Location(offset,length) );
    }
    void add_location( Location location ) { locations.push_back(location); }

    void write ( ostream & out = cout ) const;
    void clear(void);
  };


  class Sentence {
  public:
    std::map<string,string> infons;
    int offset;
    string text;
    std::vector<Annotation> annotations;
    std::vector<Relation> relations;

    void write ( ostream & out = cout ) const;
    void clear(void);
  };


  class Passage {
  public:
    std::map<string,string> infons;
    int offset;
    std::string text;                  // soon u8string;
    std::vector<Sentence> sentences;
    std::vector<Annotation> annotations;
    std::vector<Relation> relations;

    void write ( ostream & out = cout ) const;
    void clear(void);
  };


  class Document {
  public:
    string id;
    std::map<string,string> infons;
    std::vector<Passage> passages;
    std::vector<Relation> relations;

    void write ( ostream & out = cout ) const;
    void clear(void);
  };


  class Collection {
  public:
    string source;
    string date;
    string key;
    std::map<string,string> infons;
    std::vector<Document> documents;

    void write ( ostream & out = cout ) const;
    void clear(void);
  };

}

#define BIOC
#endif
