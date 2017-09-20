#include "BioC.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <sstream>

using std::cout;
using std::map;
using std::ostream;
using std::string;
using std::vector;

namespace BioC {

  void write( ostream & out, const map<string,string> & infons ) {
    map<string,string>::const_iterator iter;
    for ( iter = infons.begin(); iter != infons.end(); ++iter ) {
      out << iter->first << ": " << iter->second << '\n';
    }
  }


  ostream & operator<<( ostream & out, const Node & node ) {
    node.write(out);
    return out;
  }

  void Node::write( ostream & out ) const {
    out << refid << ':' << role;
  }


  void Relation::write ( ostream & out ) const {
    out << "id: " << id << '\n';
    BioC::write(out, infons);
    for ( int i = 0; i < nodes.size(); ++i ) 
      out << "node: " << nodes[i] << '\n';
  }

  void Relation::clear(void) {
    id.clear();
    infons.clear();
    nodes.clear();
  }

  void Location::write( ostream & out ) const {
    out << offset << ':' << length;
  }

  ostream & operator<<( ostream & out, const Location & location ) {
    location.write(out);
    return out;
  }

  int Annotation::num_segments(void) const {
    return locations.size();
  }

  int Annotation::offset( int segment ) const {
    return locations[segment].offset;
  }

  int Annotation::length( int segment ) const {
    return locations[segment].length;
  }

  void Annotation::write ( ostream & out ) const {
    out << "id: " << id << '\n';
    BioC::write(out, infons);
    for ( int i = 0; i < locations.size(); ++i )
      out << "location: " << locations[i] << '\n';
    out << "text: " << text << '\n';
  }

  void Annotation::clear(void) {
    id.clear();
    infons.clear();
    locations.clear();
    text.clear();
  }

  void  Sentence::write ( ostream & out ) const {
    BioC::write(out, infons);
    out << "offset: " << offset << '\n';
    if ( not text.empty() )
      out << "text: " << text << '\n';
    for ( int i = 0; i < annotations.size(); ++i ) 
      annotations[i].write(out);
    for ( int i = 0; i < relations.size(); ++i ) 
      relations[i].write(out);
  }

  void Sentence::clear(void) {
    infons.clear();
    offset = 0;
    text.clear();
    annotations.clear();
    relations.clear();
  }

  void Passage::write ( ostream & out ) const {
    BioC::write(out, infons);
    out << "offset: " << offset << '\n';
    if ( not text.empty() )
      out << "text: " << text << '\n';
    for ( int i = 0; i < sentences.size(); ++i ) 
      sentences[i].write(out);
    for ( int i = 0; i < annotations.size(); ++i ) 
      annotations[i].write(out);
    for ( int i = 0; i < relations.size(); ++i ) 
      relations[i].write(out);
  }

  void Passage::clear(void) {
    infons.clear();
    offset = 0;
    text.clear();
    sentences.clear();
    annotations.clear();
    relations.clear();
  }

  void Document::write ( ostream & out ) const {
    out << "id: " << id << '\n';
    BioC::write(out, infons);
    for ( int i = 0; i < passages.size(); ++i ) 
      passages[i].write(out);
    for ( int i = 0; i < relations.size(); ++i ) 
      relations[i].write(out);
  }

  void Document::clear(void) {
    id.clear();
    infons.clear();
    passages.clear();
    relations.clear();
  }

  void Collection::write ( ostream & out ) const {
    out << "source: " << source << '\n';
    out << "date: " << date << '\n';
    out << "key: " << key << '\n';
    BioC::write(out, infons);
    for ( int i = 0; i < documents.size(); ++i ) 
      documents[i].write(out);
  }

  void Collection::clear(void) {
    source.clear();
    date.clear();
    key.clear();
    infons.clear();
    documents.clear();
  }

}
