#include "BioC_util.hpp"

using namespace BioC;

void Node_Converter::convert( const Relation & relation,
                                      Relation & new_relation ) {
  new_relation = relation;
  /*
  new_relation.id = relation.id;
  new_relation.infons = relation.infons;
  new_relation.nodes = relation.nodes;
  */
}

void Node_Converter::convert( const Annotation & annotation,
                                      Annotation & new_annotation ) {
  new_annotation = annotation;
  /*
  new_annotation.id = annotation.id;
  new_annotation.infons = annotation.infons;
  new_annotation.locations = annotation.locations;
  new_annotation.text = annotation.text;
  */
}

void Node_Converter::convert( const Sentence & sentence,
                                      Sentence & new_sentence ) {
  new_sentence.infons = sentence.infons;
  new_sentence.offset = sentence.offset;
  new_sentence.text = sentence.text;

  for ( int i = 0; i < sentence.annotations.size(); ++i ) {
    Annotation new_annotation;
    convert( sentence.annotations[i], new_annotation );
    new_sentence.annotations.push_back( new_annotation );
  }

  for ( int i = 0; i < sentence.relations.size(); ++i ) {
    Relation new_relation;
    convert( sentence.relations[i], new_relation );
    new_sentence.relations.push_back( new_relation );
  }
}

void Node_Converter::convert( const Passage & passage,
                                      Passage & new_passage ) {
  new_passage.infons = passage.infons;
  new_passage.offset = passage.offset;
  new_passage.text = passage.text;

  for ( int i = 0; i < passage.sentences.size(); ++i ) {
    Sentence new_sentence;
    convert( passage.sentences[i], new_sentence );
    new_passage.sentences.push_back( new_sentence );
  }

  for ( int i = 0; i < passage.annotations.size(); ++i ) {
    Annotation new_annotation;
    convert( passage.annotations[i], new_annotation );
    new_passage.annotations.push_back( new_annotation );
  }

  for ( int i = 0; i < passage.relations.size(); ++i ) {
    Relation new_relation;
    convert( passage.relations[i], new_relation );
    new_passage.relations.push_back( new_relation );
  }
}

void Node_Converter::convert( const Document & document,
                                      Document & new_document ) {
  new_document.id = document.id;
  new_document.infons = document.infons;
    
  for ( int i = 0; i < document.passages.size(); ++i ) {
    Passage new_passage;
    convert( document.passages[i], new_passage );
    new_document.passages.push_back( new_passage );
  }

  for ( int i = 0; i < document.relations.size(); ++i ) {
    Relation new_relation;
    convert( document.relations[i], new_relation );
    new_document.relations.push_back( new_relation );
  }
}

void Node_Converter::convert( const Collection & collection,
                                      Collection & new_collection ) {
  new_collection.source = collection.source;
  new_collection.date = collection.date;
  new_collection.key = collection.key;
  new_collection.infons = collection.infons;
  
  for ( int i = 0; i < collection.documents.size(); ++i ) {
    Document new_document;
    convert( collection.documents[i], new_document );
    new_collection.documents.push_back( new_document );
  }
}
