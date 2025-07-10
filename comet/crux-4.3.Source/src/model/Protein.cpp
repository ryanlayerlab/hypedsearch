/*************************************************************************//**
 * \file protein.cpp
 * \brief Object for representing a single protein.
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <vector>
#include "util/Params.h"
#include "util/utils.h"
#include "util/crux-utils.h"
#include "parameter.h"
#include "util/GlobalParams.h"
#include "Peptide.h"
#include "Protein.h"
#include "PeptideSrc.h"
#include "Database.h"
#include "io/carp.h"
#include "PeptideConstraint.h"
#include "ProteinPeptideIterator.h"

using namespace std;
using namespace Crux;
/**
 * Constants
 */
static const int PROTEIN_ID_LENGTH = 100;
static const int PROTEIN_SEQUENCE_LENGTH = 40000;
static const int PROTEIN_ANNOTATION_LENGTH = 100;
static const int LONGEST_LINE  = PROTEIN_ID_LENGTH + PROTEIN_ID_LENGTH;
static const int FASTA_LINE = 50;
static const int SMALLEST_MASS = 57;
static const int LARGEST_MASS = 190;

/**
 * Set member variables to default values.
 */
void Protein::init() {
  database_ = NULL;
  offset_ = 0;
  protein_idx_ = 0;
  is_light_ = false;
  is_memmap_ = false;
  id_.clear();
  sequence_ = NULL;
  length_ = 0;
  annotation_.clear();
}

/**
 * \returns An (empty) protein object.
 */
Protein::Protein() {
  init();
}

/**
 * \returns A new protein object(heavy).
 * The protein is does not constain a database, users must provide one.
 */
Protein::Protein(
  const char*         id, ///< The protein sequence id. -in
  const char*   sequence, ///< The protein sequence. -in
  unsigned int length, ///< The length of the protein sequence. -in
  const char* annotation,  ///< Optional protein annotation.  -in
  unsigned long int offset, 
  ///< The file location in the source file in the database -in
  unsigned int protein_idx, ///< The index of the protein in it's database.-in  
  Database* database ///< the database of its origin
  )
{
  init();
  setId(id);
  setSequence(sequence);
  setLength(length);
  setAnnotation(annotation);
  setOffset(offset);
  setProteinIdx(protein_idx);
  setIsLight(false);
  database_ = Database::copyPtr(database); 
  is_memmap_ = false;
}         




/**
 * \returns A new light protein object.
 */
Protein* Protein::newLightProtein(
  unsigned long int offset, 
  ///< The file location in the source file in the database -in
  unsigned int protein_idx ///< The index of the protein in it's database. -in
  )
{
  Protein* protein = new Protein();
  protein->setIsLight(true);
  protein->setOffset(offset);
  protein->setProteinIdx(protein_idx);
  return protein;
}


/**
 * convert light protein to heavy, by parsing all the sequence from fasta file
 * \returns TRUE if successfully converts the protein to heavy 
 */
bool Protein::toHeavy()
{
  // protein already heavy
  if(!is_light_){
    return true;
  }
  
  FILE* file = database_->getFile();
  
  // rewind to the begining of the protein to include ">" line
  fseek(file, offset_, SEEK_SET);

  // failed to parse the protein from fasta file
  // protein offset is set in the parse_protein_fasta_file method
  if(!parseProteinFastaFile(file)){
    carp(CARP_ERROR, 
         "failed convert protein to heavy, cannot parse fasta file");
    return false;
  }
      
  is_light_ = false;
  
  return true;
}                            

/**
 * covert heavy protein back to light
 * \returns TRUE if successfully converts the protein to light
 */
bool Protein::toLight()
{
  // protein already light
  if(is_light_){
    return true;
  }
  // free all char* in protein object
  free(sequence_);
  sequence_ = NULL;
  annotation_.clear(); 
  id_.clear();
  return (is_light_ = true);
}                            

/**
 * Frees an allocated protein object.
 */
Protein::~Protein() 
{
  // FIXME what is the point of checking this?
  if (!is_memmap_ && !is_light_) {
    if (sequence_ != NULL) {
      free(sequence_);
    }
  }
}

/**
 * Prints a protein object to file.
 * if light protein coverts it to heavy protein
 */
void Protein::print(
  FILE* file ///< output stream -out
  )
{
  // covnert to heavy protein
  if (is_light_) {
    toHeavy();
  }
  int   sequence_index;
  int   sequence_length = getLength();
  char* sequence = getSequence();
  string& id = getIdPointer();
  string annotation = getAnnotation();
  
  fprintf(file, ">%s %s\n", id.c_str(), annotation.c_str());

  sequence_index = 0;
  while (sequence_length - sequence_index > FASTA_LINE) {
    fprintf(file, "%.*s\n", FASTA_LINE, &(sequence[sequence_index]));
    sequence_index += FASTA_LINE;
  }
  fprintf(file, "%s\n\n", &(sequence[sequence_index]));

  free(sequence);
}

/**
 * \returns the starting location of the sequence in a protein.  If not found, returns -1
 */
int Protein::findStart(
  string peptide_sequence, ///< the sequence to find
  string prev_aa, ///< the previous amino acid for the sequence
  string next_aa ///< the next amino acid for the sequence
) {
  if (getSequencePointer() == NULL) {
    return -1;
  }
  string protein_seq = getSequencePointer();
  if (prev_aa == "-" && StringUtils::StartsWith(protein_seq, peptide_sequence)) {
    return 1;
  } else if (next_aa == "-" && StringUtils::EndsWith(protein_seq, peptide_sequence)) {
    return getLength() - peptide_sequence.length();
  }
  // find flankN + sequence + flankC
  // prev_aa and next_aa may be "" if no flanking aa's were input.
  string seq = prev_aa + peptide_sequence + next_aa;
  size_t pos = protein_seq.find(seq);
  bool flanks_exist = prev_aa != "" && next_aa != "";
  if (pos != string::npos && flanks_exist) {
    return pos + 2;
  }
  // failed, find sequence
  seq = peptide_sequence;
  if ((pos = protein_seq.find(seq)) != string::npos) {
    return pos + 1;
  }
  // failed, find flankN + sequence + flankC with I/L equivalence
  seq = prev_aa + peptide_sequence + next_aa;
  std::replace(seq.begin(), seq.end(), 'L', 'I');
  std::replace(protein_seq.begin(), protein_seq.end(), 'L', 'I');
  if ((pos = protein_seq.find(seq)) != string::npos && flanks_exist) {
    return pos + 2;
  }
  // failed, find sequence with I/L equivalence
  seq = peptide_sequence;
  std::replace(seq.begin(), seq.end(), 'L', 'I');
  if ((pos = protein_seq.find(seq)) != string::npos) {
    return pos + 1;
  }

  // complain if the peptide is truly required
  if ( Params::GetBool("find-peptides") ) {
    carp(CARP_ERROR, "could not find %s in protein %s\n%s", peptide_sequence.c_str(), getIdPointer().c_str(), getSequencePointer());
  }
  return 0;
}

bool Protein::isPostProcess() {
  return false;
}

/**
 * prints a binary representation of the protein
 * 
 * FORMAT (no line break)
 * <int: id length><char: id><int: annotation length>
    <char: annotation><int: sequence length><char: sequence>
 *
 * Make sure when reading the binary data, add one to the length so
 * that it will read in the terminating char as well.
 */
void Protein::serialize(
  FILE* file ///< output stream -out
  )
{
  // covnert to heavy protein
  if (is_light_) {
    toHeavy();
  }
  
  int id_length = id_.length();
  int annotation_length = annotation_.length();

  // write the protein id length
  fwrite(&id_length, sizeof(int), 1, file);
  
  // write the protein id 
 // include "/0"
  fwrite(id_.c_str(), sizeof(char), id_length+1, file);

  // write the protein annotation length
  fwrite(&annotation_length, sizeof(int), 1, file);

  // write the protein annotation
  // include "/0"
  fwrite(annotation_.c_str(), sizeof(char), annotation_length+1, file);
  
  // write the protein sequence length
  fwrite(&length_, sizeof(unsigned int), 1, file);
  
  // write the protein sequence
  // include "/0"
  fwrite(sequence_, sizeof(char), length_+1, file);
}


/**
 * Copies protein object src to dest.
 * assumes that the protein is heavy
 * dest must be a heap allocated object 
 */
void Protein::copy(
  Protein* src,///< protein to copy -in
  Protein* dest ///< protein to copy to -out
  )
{
  char* sequence = src->getSequence();
  
  dest->setId(src->getIdPointer());
  dest->setSequence(sequence);
  dest->setLength(src->getLength());
  dest->setAnnotation(src->getAnnotationPointer());
  dest->setOffset(src->offset_);
  dest->setProteinIdx(src->protein_idx_);
  dest->setIsLight(src->is_light_);
  dest->database_ = src->database_;
  
  free(sequence);
}


/**
 * Parses a protein from an memory mapped binary fasta file
 * the protein_idx field of the protein must be added before or after
 * you parse the protein 
 * protein must be a heap allocated
 * 
 * Assume memmap pointer is set at beginning of protein
 * Assume protein binary format (no line break)
 * <int: id length><char: id><int: annotation length>
     <char: annotation><int: sequence length><char: sequence>
 *
 * modifies the *memmap pointer!
 * \returns TRUE if success. FALSE is failure.
 */
bool Protein::parseProteinBinaryMemmap(
  char** memmap 
  ///< a pointer to a pointer to the memory mapped binary fasta file -in
  )
{
  int id_length = 0;
  int annotation_length = 0;
  int sequence_length = 0;

  /* FIXME, maybe use this to check if still within file
  if(*memmap_as_char[0] == EOF){
    carp(CARP_ERROR, "end of file");
  }
  */

  /***set protein ID***/

  // read id length
  id_length = *((int *) *memmap);

  // reset pointer to start of id
  *memmap += sizeof(int);

  // set protein id to mem mapped id
  id_ = *memmap;

  // reset pointer to move to annotation_length
  *memmap += (id_length + 1);


  /***set protein annotation***/

  // read annotation length
  annotation_length = *((int *) *memmap);

  // reset pointer to start of annotation
  *memmap += sizeof(int);

  // set protein annotation to mem mapped annotation
  annotation_ = *memmap;

  // reset pointer to move to sequence_length
  *memmap += (annotation_length + 1);


  /***set protein sequence***/
  
  // read sequence length
  sequence_length = *((int *) *memmap);
  length_ = sequence_length;

  // reset pointer to start of sequence
  *memmap += sizeof(int);

  // set protein annotation to mem mapped sequence
  sequence_ = *memmap;

  // reset pointer to move to start of next protein
  *memmap += (sequence_length + 1);
  
  // now this protein has been created from memory mapped!
  is_memmap_ = true;

  return true;
}

// FIXME ID line and annotation might need to be fixed
VERBOSE_T verbosity = NORMAL_VERBOSE;
/**
 * Parses a protein from an open (FASTA) file.
 * the protein_idx field of the protein must be added before or after
 * you parse the protein  
 * \returns TRUE if success. FALSE is failure.
 * protein must be a heap allocated
 */
bool Protein::parseProteinFastaFile(
  FILE* file ///< fasta file -in
  )
{
  static char name[LONGEST_LINE];    ///< Just the sequence ID.
  static char desc[LONGEST_LINE];    ///< Just the comment field.
  static char buffer[PROTEIN_SEQUENCE_LENGTH];///> The sequence to read in.
  static unsigned int sequence_length; // the sequence length

  // Read the title line.
  if (!readTitleLine(file, name, desc)) {
    return(false);
  }
  
  buffer[0] = '\0';

  // Read the sequence.
  if (!readRawSequence(file, name, PROTEIN_SEQUENCE_LENGTH, buffer, &sequence_length)) {
    carp(CARP_FATAL, "Sequence %s is too long.\n", name);
    exit(1);
  }

  // update the protein object.
  setLength(sequence_length);
  setId(name);
  setSequence(buffer);
  setAnnotation(desc);

  return(true);

}

/**************************************************/

/**
 * FASTA file parsing code
 * AUTHOR: William Stafford Noble
 * modified by Chris Park
 */

/**
 * Find the beginning of the next sequence, and read the sequence ID
 * and the comment.
 */
bool Protein::readTitleLine(
  FILE* fasta_file, ///< file to read -in
  char* name, ///< write protein name here -out
  char* description) ///< write description here -out
{
  static char id_line[LONGEST_LINE];  // Line containing the ID and comment.
  int a_char;                         // The most recently read character.

  // Read until the first occurrence of ">".
  while ((a_char = getc(fasta_file)) != '>') {
    // If we hit the end of the file, return FALSE.
    if (a_char == EOF) {
      return(false);
    }  
  }
  // set protein offset                   FIXME: might not need to "-1" -CHRIS
  offset_ = ftell(fasta_file) - 1;

  /**
   * chris edited, added this block to make sure all of comment line
   * is read although might not be stored, to ensure the file* is at
   * start of the sequence
   */

  {
    char* new_line = NULL;
    int line_length;
    size_t buf_length = 0;

    if ((line_length = getline(&new_line, &buf_length, fasta_file)) == -1) {
      carp(CARP_FATAL, "Error reading Fasta file.\n");
    }
    strncpy(id_line, new_line, LONGEST_LINE-1);
    free(new_line);
  }

  // Remove EOL.
  id_line[strlen(id_line) - 1] = '\0';

  // Extract the ID from the beginning of the line.
  if (sscanf(id_line, "%s", name) != 1) {
    carp(CARP_FATAL, "Error reading sequence ID.\n%s\n", id_line);
  }

  // Store the rest of the line as the comment.
  strcpy(description, &(id_line[strlen(name)+1]));

  return(true);
}


/****************************************************************************
 * Read raw sequence until a '>' is encountered or too many letters
 * are read.  The new sequence is appended to the end of the given
 * sequence.
 *
 * Return: Was the sequence read completely?
 ****************************************************************************/
bool Protein::readRawSequence(
  FILE* fasta_file,   // Input Fasta file.
  char* name,         // Sequence ID (used in error messages).
  unsigned int   max_chars,    // Maximum number of allowed characters.
  char* raw_sequence, // Pre-allocated sequence.
  unsigned int* sequence_length // the sequence length -chris added
  )
{
  int a_char;
  unsigned int i_seq;
  bool return_value = true;

  // Start at the end of the given sequence.
  i_seq = strlen(raw_sequence);
  assert((unsigned int)strlen(raw_sequence) < max_chars);

  // Read character by character.
  while ((a_char = getc(fasta_file)) != EOF) {

    // Check for the beginning of the next sequence.
    if (a_char == '>') {
      // Put the ">" back onto the stream for the next call to find.
      ungetc(a_char, fasta_file);
      break;
    }

    // Skip non-alphabetic characters.
    if (!isalpha((int)a_char)) {
      if ((a_char != ' ') && (a_char != '\t') && (a_char != '\n') && (a_char != '\r')) {
        carp(CARP_WARNING,"Skipping character %c in sequence %s.",
             a_char, name);
      }

    } else {

      // Convert invalid characters to X.
      a_char = toupper((int)a_char);

      /**
       * Check the ASCII code.  If the char is above or below the
       * A(65)~Z(90)range, convert the character to an 'X'.
       */
      if ( (int)a_char < 65 || (int)a_char  > 90 ) {
        carp(CARP_WARNING, "Converting illegal character %c to X ",
             a_char);
        carp(CARP_WARNING, "in sequence %s.", name);
        a_char = 'X';
      }
      
      raw_sequence[i_seq] = a_char;
      i_seq++;
    }
    if (i_seq >= max_chars) {
      return_value = false;
      break;
    }
  }
  raw_sequence[i_seq] = '\0';
  *sequence_length = i_seq; // chris added

  return(return_value);
}


/**
 * end of FASTA parsing
 * Thanks Bill!
 */

/**
 * Change the sequence of a protein to be a randomized version of
 * itself.  The method of randomization is dependent on the
 * decoy_type (shuffle or reverse).  The name of the protein is also
 * changed by prefixing with "decoy-prefix"
 */
void Protein::shuffle(
  DECOY_TYPE_T decoy_type){ ///< method for shuffling
  char* decoy_str = decoy_type_to_string(decoy_type);
  carp(CARP_DEBUG, "Shuffling protein %s as %s", 
       id_.c_str(), decoy_str);
  free(decoy_str);

  switch(decoy_type){
  case NO_DECOYS:
    return;

  case PROTEIN_SHUFFLE_DECOYS:
    carp(CARP_DEBUG, "shuffling");
    shuffle_array(sequence_, strlen(sequence_));
    break;

  case PEPTIDE_SHUFFLE_DECOYS:
    peptideShuffleSequence();
    break;
    
  case INVALID_DECOY_TYPE:
  case NUMBER_DECOY_TYPES:
    carp(CARP_FATAL, "Illegal decoy type for shuffling protein.");
    break;
  }

  // change the protein name
  string prefix = Params::GetString("decoy-prefix");
  id_ = prefix + id_;
}

/** 
 * Access routines of the form get_<object>_<field> and set_<object>_<field>. 
 */

/**
 * Additional get and set methods
 */

/**
 *\returns the id of the protein
 * returns a heap allocated new copy of the id
 * user must free the return id
 * assumes that the protein is heavy
 */
string Protein::getId()
{

  if(is_light_){
    carp(CARP_FATAL, "Cannot get ID from light protein.");
  }
  
  return id_;

}

/**
 *\returns a pointer to the id of the protein
 * assumes that the protein is heavy
 */
string& Protein::getIdPointer()
{
  if(is_light_){
    carp(CARP_FATAL, "Cannot get ID pointer from light protein.");
  }
  return id_; 
}

/**
 * sets the id of the protein
 */
void Protein::setId(
  const string& id ///< the sequence to add -in
  )
{
  id_ = id;
}

/**
 *\returns the sequence of the protein
 * returns a heap allocated new copy of the sequence
 * user must free the return sequence 
 * assumes that the protein is heavy
 */
char* Protein::getSequence(
  int offset
) {
  if(is_light_){
    carp(CARP_FATAL, "Cannot get sequence from light protein.");
  }
  unsigned int sequence_length = strlen(sequence_) +1-offset; // +\0
  char * copy_sequence = 
    (char *)mymalloc(sizeof(char)*sequence_length);
  return strncpy(copy_sequence, sequence_+offset, sequence_length);  
}

/**
 *\returns a pointer to the sequence of the protein
 * assumes that the protein is heavy
 */
char* Protein::getSequencePointer(
  int offset
)
{
  if (is_light_) {
    carp(CARP_FATAL, "Cannot get sequence pointer from light protein.");
  }
  return sequence_+offset;
}

/**
 * sets the sequence of the protein
 */
void Protein::setSequence(
  const char* sequence ///< the sequence to add -in
  )
{

  free(sequence_);
  unsigned int sequence_length = strlen(sequence) +1; // +\0
  char * copy_sequence = 
    (char *)mymalloc(sizeof(char)*sequence_length);
  sequence_ =
    strncpy(copy_sequence, sequence, sequence_length);  
}

/**
 *\returns the length of the protein
 * assumes that the protein is heavy
 */
unsigned int Protein::getLength()
{
  return length_;
}

/**
 * sets the id of the protein
 */
void Protein::setLength(
  unsigned int length ///< the length to add -in
  )
{
  length_ = length;
}

/**
 *\returns the annotation of the protein
 * returns a heap allocated new copy of the annotation
 * user must free the return annotation
 * assumes that the protein is heavy
 */
string Protein::getAnnotation()
{
  if(is_light_){
    carp(CARP_FATAL, "Cannot get annotation from light protein.");
  }

  return annotation_;
}

/**
 *\returns A const pointer to the annotation of the protein.
 */
const string& Protein::getAnnotationPointer(){
  return annotation_;
}
 
/**
 * sets the annotation of the protein
 */
void Protein::setAnnotation(
  const string& annotation ///< the sequence to add -in
  )
{
  annotation_ = annotation;
}

/**
 * sets the offset of the protein in the fasta file
 */
void Protein::setOffset(
  unsigned long int offset 
  ///< The file location in the source file in the database -in
  )
{
  offset_ = offset;
}

/**
 *\returns the offset the protein
 */
unsigned long int Protein::getOffset()
{
  return offset_;
}

/**
 * sets the protein_idx (if, idx=n, nth protein in the fasta file)
 */
void Protein::setProteinIdx(
  unsigned int protein_idx ///< The index of the protein in it's database. -in
  )
{
  // carp(CARP_DETAILED_DEBUG, "set protein idx = %i", protein_idx);
  protein_idx_ = protein_idx;
}

/**
 *\returns the protein_idx field
 */
unsigned int Protein::getProteinIdx()
{
  return protein_idx_;
}

/**
 * sets the is_light field (is the protein a light protein?)
 */
void Protein::setIsLight(
  bool is_light ///< is the protein a light protein? -in
  )
{
  is_light_ = is_light;
}

/**
 *\returns TRUE if the protein is light protein
 */
bool Protein::getIsLight()
{
  return is_light_;
}

/**
 * sets the database for protein
 */
void Protein::setDatabase(
  Database*  database ///< Which database is this protein part of -in
  )
{
  database_ = Database::copyPtr(database);
}

/**
 *\returns Which database is this protein part of
 */
Database* Protein::getDatabase()
{
  return database_;
}

/** 
 * Comparison function for sorting proteins by protein id.
 */
bool Crux::protein_id_less_than(Protein* protein_one, Protein* protein_two){
  int compare = protein_one->getIdPointer().compare(protein_two->getIdPointer());
  return (compare > 0);
}

/**
 * Rearrange the sequence_ between cleavage sites, keeping residues
 * on either side of a cleavage in place.  Get enzyme from
 * parameters.  Cases of NO_ENZYME or NON_SPECIFIC_DIGEST are the same
 * as shuffling the whole protein.  Same behavior for full and partial
 * digest, min/max length/mass and missed cleavages, i.e. shuffle
 * between every cleavage site.
 */
void Protein::peptideShuffleSequence(){
  if (sequence_ == NULL ) {
    carp(CARP_WARNING, "Cannot shuffle a NULL sequence");
    return;
  }
  if (length_ < 2 ) {
    return;
  }
 
  // get the digest rule
  ENZYME_T enzyme = GlobalParams::getEnzyme();
  // cases where peptide-shuffle is really protein shuffle
  if (enzyme == NO_ENZYME
      || GlobalParams::getDigestion() == NON_SPECIFIC_DIGEST){
    this->shuffle(PROTEIN_SHUFFLE_DECOYS);
    return;
  }

  // store valid cleavage locations
  vector<int> cleave_after_here;

  // traverse the sequence to penultimate residue
  for (size_t seq_offset = 0; seq_offset < length_ - 1; seq_offset++) {
    // mark each valid location (mark x for cleavage between x and y)
    if (ProteinPeptideIterator::validCleavagePosition(sequence_ + seq_offset,
                                                      enzyme) ){
      cleave_after_here.push_back(seq_offset);
    }
  }

  // shuffle between each cleavage site (leave a and b in place)
  int start = 0; // shuffle between prot first residue and first site
  int end = -1;
  vector<int>::iterator next_position = cleave_after_here.begin();
  for (; next_position != cleave_after_here.end(); ++next_position) {
    end = *next_position;
    shuffleRegion(start, end);
    start = end + 1; // hold in place both sides of the cleavage site
  }

  // shuffle end of sequence
  end = length_ - 1;
  if (start < end ) {
    shuffleRegion(start, end);
  }
}

/**
 * Shuffle the region of the sequence between start and end, leaving
 * start and end residues in place.  Repeat up to three times if the
 * shuffled sequence doesn't change.
 */
void Protein::shuffleRegion(
  int start, ///< index of peptide start
  int end)
{///< index of last residue in peptide
  int sub_seq_length = end - start - 1;
  char* buf = new char[sub_seq_length + 1];
  if (sub_seq_length > 1) {
    carp(CARP_DETAILED_DEBUG, "Shuffle from %d to %d.", start+1, end);

    // store the sequence before shuffling not including the unmoved residues
    strncpy(buf, sequence_ + start + 1, sub_seq_length);
    buf[sub_seq_length] = '\0';

    shuffle_array(sequence_ + start + 1, sub_seq_length);

    // check to see if it changed
    bool has_changed = strncmp(buf, sequence_ + start, sub_seq_length + 2);
    // try reshuffling up to three more times
    int count = 0;
    while ( (count < 3) && (has_changed == false)) {
      shuffle_array(sequence_ + start + 1, sub_seq_length);
      has_changed = strncmp(buf, sequence_ + start, sub_seq_length + 2);
      count++;
    }
    if (!has_changed) {
        carp(CARP_WARNING, "Unable to generate a shuffled sequence "
             "different than the original for sequence %s of protein %s "
             "at position %d.", buf, id_.c_str(), start);
    }
  }
  delete [] buf;
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

