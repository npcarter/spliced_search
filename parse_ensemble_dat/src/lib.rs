// parse_ensemble_dat: takes an ensemble dat file that describes a genome and
// parses it into DNA chunks that contain proteins along with vectors that
// describe their introns, exons, start points, and end points.
use std::fs::File;
use std::iter::Peekable;
use std::io::{BufReader, Lines};
use std::cmp::{min, max};
use rand::{RngExt};
// For exons, need to know where they start and end,
// plus whether they are complementary to the DNA sequence.
// start/end are initially absolute positions within the contig, 
// but get converted to offsets within the DNA region of the protein
pub struct Exon{
    pub start: i64,
    pub end: i64,
    pub complement: bool,
}

// For introns, need to know where they start and end in the contig.
pub struct Intron{
    pub start: i64,
    pub end: i64,
}


pub struct Protein{
    pub name: String,
    pub uniprot_name: String,
    pub mrna_start: i64, // Start of the mRNA within the contig
    pub mrna_end: i64, // End of the mRNA within the contig
    pub coding_start: i64, //start position of the coding sequence within the mRNA
    pub coding_end: i64, // end position of the coding sequence within the mRNA
    pub introns: Vec<Intron>,
    pub exons: Vec<Exon>,
    pub codon_start: i8, // How many positions (0-2) is the first codon start offset from the
    // coding sequence start
    // use these because we want to pad the start and end of the transcript for learning.
    pub dna_start: i64, // start position of the DNA sequence within the contig
    pub dna_end: i64, // end position of the DNA sequence within the contig
    pub dna_sequence: Vec<char>,  //DNA sequence of the protein
    pub protein_sequence: Vec<char>,
    pub translation_table: i8,
}

pub struct Gene{
    pub name: String,
  //  dna_sequence: Vec<char>,
    pub proteins: Vec<Protein>,
}

pub struct Contig{
    pub length: i64,
    pub name: String,
    pub organism: String,
    pub taxonomy: Vec<String>,
    pub dna_sequence: Vec<char>,
    pub genes: Vec<Gene>,
}

// find the complement of a residue
pub fn complement(residue: char) -> char{
    match residue {
        'A' => 'T',
        'C' => 'G',
        'G' => 'C',
        'T' => 'A',
        'N' => 'N',
        'R' => 'R',
        'Y' => 'Y',
        'K' => 'K',
        'M' => 'M',
        'S' => 'S',
        'W' => 'W',
        'B' => 'B',
        'D' => 'D',
        'H' => 'H',
        'V' => 'V',
        _ => panic!("Invalid residue {}", residue),
    }
}

pub struct Gap{
    pub start: i64,
    pub end: i64,
    pub name: String,
}

// translate a codon into an amino acid
pub fn translate_codon(codon: &Vec<char>) -> char{
    if codon.len() != 3 { panic!("Codon must be three characters long"); }

    let amino_acid:char;
    let codon_string = codon.iter().collect::<String>();
    match codon_string.as_str(){
        "AAA" => amino_acid = 'K',
        "AAC" => amino_acid = 'N',
        "AAG" => amino_acid = 'K',
        "AAT" => amino_acid = 'N',
        "ACA" => amino_acid = 'T',
        "ACC" => amino_acid = 'T',
        "ACG" => amino_acid = 'T',
        "ACT" => amino_acid = 'T',
        "AGA" => amino_acid = 'R',
        "AGC" => amino_acid = 'S',
        "AGG" => amino_acid = 'R',
        "AGT" => amino_acid = 'S',
        "ATA" => amino_acid = 'I',
        "ATC" => amino_acid = 'I',
        "ATG" => amino_acid = 'M',
        "ATT" => amino_acid = 'I',
        "CAA" => amino_acid = 'Q',
        "CAC" => amino_acid = 'H',
        "CAG" => amino_acid = 'Q',
        "CAT" => amino_acid = 'H',
        "CCA" => amino_acid = 'P',
        "CCC" => amino_acid = 'P',
        "CCG" => amino_acid = 'P',
        "CCT" => amino_acid = 'P',
        "CGA" => amino_acid = 'R',
        "CGC" => amino_acid = 'R',
        "CGG" => amino_acid = 'R',
        "CGT" => amino_acid = 'R',
        "CTA" => amino_acid = 'L',
        "CTC" => amino_acid = 'L',
        "CTG" => amino_acid = 'L',
        "CTT" => amino_acid = 'L',
        "GAA" => amino_acid = 'E',
        "GAC" => amino_acid = 'D',
        "GAG" => amino_acid = 'E',
        "GAT" => amino_acid = 'D',
        "GCA" => amino_acid = 'A',
        "GCC" => amino_acid = 'A',
        "GCG" => amino_acid = 'A',
        "GCT" => amino_acid = 'A',
        "GGA" => amino_acid = 'G',
        "GGC" => amino_acid = 'G',
        "GGG" => amino_acid = 'G',
        "GGT" => amino_acid = 'G',
        "GTA" => amino_acid = 'V',
        "GTC" => amino_acid = 'V',
        "GTG" => amino_acid = 'V',
        "GTT" => amino_acid = 'V',
        "TAA" => amino_acid = '*',
        "TAC" => amino_acid = 'Y',
        "TAG" => amino_acid = '*',
        "TAT" => amino_acid = 'Y',
        "TCA" => amino_acid = 'S',
        "TCC" => amino_acid = 'S',
        "TCG" => amino_acid = 'S',
        "TCT" => amino_acid = 'S',
        "TGA" => amino_acid = '*',
        "TGC" => amino_acid = 'C',
        "TGG" => amino_acid = 'W',
        "TGT" => amino_acid = 'C',
        "TTA" => amino_acid = 'L',
        "TTC" => amino_acid = 'F',
        "TTG" => amino_acid = 'L',
        "TTT" => amino_acid = 'F',
        _ => amino_acid = 'X',
    }

    amino_acid
}

// check whether a codon matches one of the possible start codons
pub fn start_codon(codon: &Vec<char>) -> bool{
    if codon.len() != 3 { panic!("Codon must be three characters long"); }
    let codon_string = codon.iter().collect::<String>();
    match codon_string.as_str(){
        "ATG" => true,
        "TTG" => true,
        "CTG" => true,
        _ => false,
    }
}

// translate a DNA sequence into a protein sequence
pub fn translate_sequence(mrna_sequence: &Vec<char>, codon_start: i8) -> Vec<char>{
    let mut protein_sequence:Vec<char> = Vec::new();
    let dna_sequence = mrna_sequence[codon_start as usize ..].to_vec();
   // if dna_sequence.len() % 3 != 0 { panic!("DNA sequence must be a multiple of three characters long"); }
   // special-case the start codon
   if start_codon(&dna_sequence[0..3].to_vec()){
    protein_sequence.push('M');
   }
   else{
    protein_sequence.push(translate_codon(&dna_sequence[0..3].to_vec()));
   }
    for i in 1..(dna_sequence.len()/3){
        protein_sequence.push(translate_codon(&dna_sequence[i*3..(i*3)+3].to_vec()));
    }

    protein_sequence
}




// Because the dat file gives the DNA sequence for the contig at the end of its
// description, after we've parsed the contig, we need to go through and fill in
// the protein sequences for each gene.  This function handles that.
pub fn parse_contig(contig: &mut Contig){
  //  println!("Contig {} has length {} and organism {}", contig.name, contig.length, contig.organism);
    let mut rng = rand::rng(); // 
    for gene in contig.genes.iter_mut() {
        for protein in gene.proteins.iter_mut() {
            let mut check_sequence: Vec<char> = Vec::new();
            // pad the start and end of the transcript to a random length between 100 and 10% of the transcript length
            let pad_max_length = max((protein.mrna_end - protein.mrna_start + 1) /10, 101);
            // don't pad more than the start or end of the transcript
            // subtract one from mrna_start because the positions are 1-indexed and we don't want to get a dna_start = 0
            // if the mrna_start is near the beginning of the contig
            let padding_start = min(rng.random_range(100..pad_max_length), protein.mrna_start-1);

            // don't subtract one here because we want to pad the mrna by zero positions if it ends at the end of the contig
            let padding_end = min(rng.random_range(100..pad_max_length), contig.dna_sequence.len() as i64 - protein.mrna_end);
            protein.dna_start = protein.mrna_start - padding_start;
            protein.dna_end = protein.mrna_end + padding_end;
            if !protein.exons[0].complement {  // forward strand sequence
                protein.coding_start += padding_start; // adjust the location of the coding region to account for the padding
                protein.coding_end += padding_start;
                for index in protein.dna_start..protein.dna_end + 1
                { protein.dna_sequence.push(contig.dna_sequence[index as usize - 1]); }

                /*     //need to handle first exon separately because coding sequence may not start at start of exon
                     for index in (protein.coding_start + protein.exons[0].start)..protein.exons[0].end+1{
                         check_sequence.push(protein.dna_sequence[(index -protein.mrna_start) as usize ]);
                     }
     */
                for exon in protein.exons.iter_mut() {
                    assert!(exon.start >= protein.dna_start); // sanity-check the exon start/end positions
                    assert!(exon.end >= protein.dna_start);
                    exon.start -= protein.dna_start; // Convert from index within contig to offset within the protein's DNA region
                    exon.end -= protein.dna_start;
                    for index in max(protein.coding_start,exon.start)..(min(exon.end, protein.coding_end) + 1) {
                        check_sequence.push(protein.dna_sequence[index as usize]);
                    }
                    //println!("{}, {}, length={}", exon.start-protein.mrna_start, min(exon.end, protein.coding_end+protein.mrna_start)-protein.mrna_start, check_sequence.len());
                }
            } else { //reverse strand sequence
                //compute the offsets to the coding start and end from the end of the sequence to handle reversing, then swap them 
                //because the end of the unreversed seqence is the start of the reversed and vice-versa
                let reverse_end_offset = (protein.dna_end - protein.mrna_start) - protein.coding_start;  // offset in reversed protein = DNA length - offset from beginning
                let reverse_start_offset = (protein.dna_end - protein.mrna_start) - protein.coding_end;
                protein.coding_start = reverse_start_offset;
                protein.coding_end = reverse_end_offset;
                for index in (protein.dna_start..protein.dna_end + 1).rev()
                {
                    protein.dna_sequence.push(complement(contig.dna_sequence[index as usize - 1]));
                }
                //println!("{}: ", protein.name);
                for exon in protein.exons.iter_mut() {
                    assert!(exon.start <= protein.dna_end); // sanity-check the exon start/end positions
                    assert!(exon.end <= protein.dna_end);
                    exon.start = protein.dna_end - exon.start; // Convert from index within contig to offset within the protein's DNA region
                    exon.end = protein.dna_end - exon.end; // accounting for the fact that we reverse direction of the strand
 //                   println!("start: {}, end: {}", max(exon.end, protein.dna_end - (protein.coding_end +protein.dna_start)), min(protein.dna_end - (protein.coding_start+protein.dna_start), exon.start)+1);

                    // need to convert coding_start, coding_end back into absolute indices within transcript to get delta from protein.dna_end
                    for index in max(exon.end, protein.coding_start)..min(protein.coding_end, exon.start)+ 1 {
                        // We've already reversed the protein sequence, so exon.end (as an offset into protein.dna_sequence) will be less than exon.start
                        check_sequence.push(protein.dna_sequence[(index) as usize]);
                       //print!("{} ", protein.mrna_end - index);
                       //println!("{} {}",(protein.mrna_end - index),protein.dna_sequence[(protein.mrna_end- index)as usize]);
                    }
                }
            }
            if protein.translation_table ==0{ // don't support alternate translations yet, but this is just self-checks
                let mut translated_sequence = translate_sequence(&check_sequence, protein.codon_start);
                if translated_sequence.last().unwrap() == &'*' { translated_sequence.pop(); } // handle trailing stop codon
                let mut found_missmatch = false;

                let mut start_offset = 0;
                if protein.protein_sequence.len() > translated_sequence.len(){ // see if we need to skip leading 'X' characters in the ENSEMBL sequence
                    while protein.protein_sequence[start_offset] == 'X'{
                        start_offset += 1;
                    }
                }
            
                if translated_sequence.len() + start_offset != protein.protein_sequence.len() {
                    found_missmatch = true;
                }
                for i in 0..translated_sequence.len() {
                    if translated_sequence[i] != protein.protein_sequence[i+start_offset] && translated_sequence[i] != 'X' && protein.protein_sequence[i+start_offset] != 'X' && !(translated_sequence[i] == '*' && protein.protein_sequence[i+start_offset] == 'U'){
                        found_missmatch = true;
                        println!("missmatch at position {}: {} vs {}", i, translated_sequence[i], protein.protein_sequence[i]);
                        if i == 0{
                            println!("possible non-canonical start codon of {}{}{}", check_sequence[0], check_sequence[1], check_sequence[2]);
                        }
                    }
                }
                if found_missmatch{
                    println!("Error: protein {} sequence {} does not match translated sequence {}", protein.name, protein.protein_sequence.iter().collect::<String>(), translated_sequence.iter().collect::<String>());
                    let first_mismatch = translated_sequence[1..].iter()
                        .zip(protein.protein_sequence[1..].iter())
                        .position(|(a, b)| a != b);
                    if first_mismatch != None{
                        println!("First mismatch at position {}.", first_mismatch.unwrap());
                    }
                    else if protein.protein_sequence.len() > translated_sequence.len(){
                        println!("Original sequence was longer than translated sequence, {} vs {}", protein.protein_sequence.len(), translated_sequence.len());
                    }
                    else{
                        println!("unknown error in checking sequence translation");
                    }
                }
            }
        }
    }
}

// Parse the DNA sequence from the dat file
pub fn parse_sequence(fields: Vec<&str>, lines: &mut Peekable<Lines<BufReader<File>>>) -> Vec<char>{
    let mut sequence:Vec<char> = Vec::new();
    let length:i64 = fields[2].parse::<i64>().unwrap();
    let mut residue_count:i64 = 0;
    let mut next_line = match lines.peek(){
        Some(Ok(line)) => line.to_string(),
        _ => panic!("Unexpected end of file")
    };
    while next_line != "//".to_string(){
        let line = (&mut *lines).next().unwrap().unwrap();
        let fields = line.split_whitespace().collect::<Vec<&str>>();
        let end:i64 = fields.last().unwrap().parse::<i64>().unwrap();
        for residues in fields[0..(fields.len()-1)].iter(){
            for residue in residues.chars() {
                sequence.push(residue);
                residue_count += 1;
            }
        }
        if residue_count != end{println!("Error: residue count {} does not match end position {}", residue_count, end);}
        next_line = match lines.peek(){
            Some(Ok(line)) => line.to_string(),
            _ => panic!("Unexpected end of file")
        };
    }
 //   println!("Sequence length is {}", sequence.len());
    assert_eq!(length, sequence.len() as i64, "length of sequence did not match parsed sequence length");
    sequence
}

// Parse a gene from the dat file.
// a gene can contain multiple proteins, often because of alternative splicing
pub fn parse_gene(lines: &mut Peekable<Lines<BufReader<File>>>) -> Gene{
    //first, grab the start and end positions of the gene

    //let start:i64 = start_end[0].parse::<i64>().unwrap();
    //let end:i64 = start_end[1].parse::<i64>().unwrap();

    let mut the_gene = Gene{proteins: Vec::new(), name: "".to_string()};

    // now go over the remaining lines and create one or more protein descriptions
    let mut next_line = match lines.peek(){
        Some(Ok(line)) => line.to_string(),
        _ => panic!("Unexpected end of file")
    };
    let mut found_mrna = false; // gene must have at least one mRNA, track this
    // to avoid gene descriptions with "gene" in the name
    let mut next_fields = next_line.split_whitespace().collect::<Vec<&str>>();
    while next_fields[0] != "XX" && next_fields[0] != "SQ" && next_fields[0] != "//" &&
        !(next_fields[0] == "FT" && next_fields[1] == "gene" && found_mrna){
        // We're still in the gene

        let line = (&mut *lines).next().unwrap().unwrap();
        let line_fields = line.split_whitespace().collect::<Vec<&str>>();
        if line_fields[1].starts_with("/gene") { the_gene.name = line_fields[1].trim_start_matches("/gene=").trim_start_matches('\"').to_string().trim_end_matches('\"').to_string();
        }
        else {
            if line_fields[0] == "FT" && line_fields[1] == "misc_RNA" {
                found_mrna = true; // some genes just have misc RNAs and no coding regions
            }
            if line_fields[0] == "FT" && line_fields[1] == "mRNA" {
                found_mrna = true;
                the_gene.proteins.push(parse_protein(line_fields, lines));
            }
        }

        next_line = match lines.peek(){
            Some(Ok(line)) => line.to_string(),
            _ => panic!("Unexpected end of file")
        };
        next_fields = next_line.split_whitespace().collect::<Vec<&str>>();
    }
    the_gene
}

// parse the description of a protein from the dat file
pub fn parse_protein(fields: Vec<&str>, lines: &mut Peekable<Lines<BufReader<File>>>) -> Protein{
    let mut the_protein = Protein{name: "".to_string(), uniprot_name: "".to_string(), mrna_start: 0, mrna_end: 0, translation_table: 0, coding_start: 0, coding_end: 0,
        introns: Vec::new(), exons: Vec::new(), codon_start: 0, dna_start: 0, dna_end: 0, dna_sequence: Vec::new(), protein_sequence: Vec::new()};

    // First, construct the list of exons from the join statement
    for exon_string in fields[2].trim_start_matches("join(").trim_end_matches(")").trim_end_matches(',').split(',') {
        the_protein.exons.push(parse_exon(exon_string.trim_start_matches('(').trim_end_matches(')')));
    }
    let mut next_line = match lines.peek(){
        Some(Ok(line)) => line.to_string(),
        _ => panic!("Unexpected end of file")
    };
    let mut next_fields = next_line.split_whitespace().collect::<Vec<&str>>();
    while !(next_fields[0] == "FT" && next_fields[1].starts_with("/gene")){
        let line = (&mut *lines).next().unwrap().unwrap();
        // Exon start/end positions start out as absolute indices within the contig.  Later, we'll transform them into offsets within the region of DNA
        // that contains the protein.  We do this as a two-step process to handle the padding that we add to the start/end of the DNA later.
        for exon_string in line.split_whitespace().skip(1).collect::<Vec<&str>>().join("").trim_end_matches(")").trim_end_matches(',').split(',') {
            the_protein.exons.push(parse_exon(exon_string.trim_start_matches('(').trim_end_matches(')')));
        }
        next_line = match lines.peek(){
            Some(Ok(line)) => line.to_string(),
            _ => panic!("Unexpected end of file")
        };
        next_fields = next_line.split_whitespace().collect::<Vec<&str>>();
    }
    assert!(the_protein.exons.len() > 0, "Protein had no exons"); // need to have at least one exon
    // get the start and end positions of the mRNA from the exons
    if the_protein.exons[0].complement { // complement exons listed in reverse order
        the_protein.mrna_start = the_protein.exons.last().unwrap().start;
        the_protein.mrna_end = the_protein.exons.first().unwrap().end;
    }
    else{
        the_protein.mrna_start = the_protein.exons.first().unwrap().start;
        the_protein.mrna_end = the_protein.exons.last().unwrap().end;
    }

    // Now that we have the list of exons, construct the list of introns from the gaps between exons
    for i in 1..the_protein.exons.len(){
        assert_eq!(the_protein.exons[i].complement, the_protein.exons[i-1].complement);  // all exons in a protein should be either complementary or non-complementary
        if !the_protein.exons[i].complement{
        the_protein.introns.push(Intron{start: the_protein.exons[i-1].end+1, end: the_protein.exons[i].start-1});}
        else{
            the_protein.introns.push(Intron{start: the_protein.exons[i].end+1, end: the_protein.exons[i-1].start-1});
        }
        assert!(the_protein.introns.last().unwrap().end >= the_protein.introns.last().unwrap().start);

    }
    assert_eq!(the_protein.introns.len(), the_protein.exons.len()-1);

    let _skip_line = (&mut *lines).next().unwrap().unwrap(); // skip the /gene line

    let name_line = (&mut *lines).next().unwrap().unwrap();  // get the /standard_nme line
    let name_fields = name_line.split_whitespace().collect::<Vec<&str>>();
    assert!(name_fields[0] == "FT" && name_fields[1].starts_with("/standard_name"));
    the_protein.name = name_fields[1].split('"').collect::<Vec<&str>>()[1].to_string();

    // parse the CDS data to get the protein start and end positions within the gene
    let mut cds_line = (&mut *lines).next().unwrap().unwrap();
    let mut cds_fields = cds_line.split_whitespace().collect::<Vec<&str>>();
    assert!(cds_fields[0] == "FT" && cds_fields[1].starts_with("CDS"), "Misformatted CDS line: {}", cds_line);
    cds_fields.drain(0..2); // drop the FT and CDS fields
    let cds_string = cds_fields.join("");
    let mut cds_exons = cds_string.trim_start_matches("join(").trim_end_matches(")").trim_end_matches(',').split(',').collect::<Vec<&str>>();
    let first_cds_exon = parse_exon(cds_exons[0].trim_start_matches('(').trim_end_matches(')'));
    let mut last_cds_exon = parse_exon(cds_exons[cds_exons.len()-1].trim_start_matches('(').trim_end_matches(')'));

    // now, see if there are any other lines of exons
    cds_line = (&mut *lines).next().unwrap().unwrap();
    cds_fields = cds_line.split_whitespace().collect::<Vec<&str>>();
    while !cds_fields[1].starts_with("/"){
        cds_fields.drain(0..1);
        let cds_string2 = cds_fields.join("");
        cds_exons = cds_string2.trim_end_matches(")").trim_end_matches(',').split(',').collect();
   //     cds_fields[1].trim_start_matches("join(").trim_end_matches(")").trim_end_matches(',').split(',').collect::<Vec<&str>>();
        last_cds_exon = parse_exon(cds_exons[cds_exons.len()-1].trim_start_matches('(').trim_end_matches(')'));
        cds_line = (&mut *lines).next().unwrap().unwrap();
        cds_fields = cds_line.split_whitespace().collect::<Vec<&str>>();
    }
    if cds_fields[1].starts_with("/codon_start="){
        //codons are offset from the first residue in the coding region
        the_protein.codon_start = cds_fields[1].trim_start_matches("/codon_start=\"").trim_end_matches("\"").parse::<i8>().expect("Failed to parse codon_start field") -1;
        // -1 because codon_start of 1 means zero offset from the start of the sequence
    }
    else if cds_fields[1].starts_with("/transl_table="){
        the_protein.translation_table = cds_fields[1].trim_start_matches("/transl_table=").parse::<i8>().expect("Failed to parse transl_table field");
    }
    if !first_cds_exon.complement{
        the_protein.coding_start = first_cds_exon.start - the_protein.mrna_start;
        the_protein.coding_end = last_cds_exon.end - the_protein.mrna_start;
    }
    else{
        the_protein.coding_start =  last_cds_exon.start- the_protein.mrna_start;
        the_protein.coding_end = first_cds_exon.end - the_protein.mrna_start;
    }
    //println!("Protein {} has mRNA from {} to {} and coding sequence from {} to {}", the_protein.name,
      //     the_protein.mrna_start, the_protein.mrna_end, the_protein.coding_start, the_protein.coding_end);

    // parse the rest of the mRNA sequence
    let mut in_translation = false;
    // go through the rest of the protein description and find its name and translation
    while !(next_fields[0] == "FT" && (next_fields[1] == "gene" || next_fields[1] == "mRNA")) && next_fields[0] != "//"
    && next_fields[0] != "XX" && next_fields[0] != "SQ" && !(next_fields[0] == "FT" && next_fields[1] == "exon"){
        // still in this protein
        let line = (&mut *lines).next().unwrap().unwrap();
        let line_fields = line.split_whitespace().collect::<Vec<&str>>();
        if in_translation { // handle multi-line translations
            the_protein.protein_sequence.extend(line_fields[1].trim_end_matches('"').chars());
        }
        if line_fields[0] == "FT" && line_fields[1].starts_with("/translation"){ the_protein.protein_sequence = line_fields[1].split('"').collect::<Vec<&str>>()[1].chars().collect::<Vec<char>>();
        in_translation = true;}
        if line_fields[0] == "FT" && line_fields[1].starts_with("/db_xref=\"Uniprot/"){ the_protein.uniprot_name = line_fields[1].split('/').collect::<Vec<&str>>()[2].to_string().trim_end_matches('"').to_string();}

        next_line = match lines.peek(){
            Some(Ok(line)) => line.to_string(),
            _ => panic!("Unexpected end of file")
        };
        next_fields = next_line.split_whitespace().collect::<Vec<&str>>();
    }
    //println!("Protein {} has name {} and translation {}", the_protein.name, the_protein.uniprot_name, the_protein.protein_sequence.iter().collect::<String>());
    the_protein
}


// generate an exon from the string that describes it.
pub fn parse_exon(exon_string: &str) -> Exon{
    let start_end:Vec<&str>;
    let complement:bool = exon_string.starts_with("complement");
    if !complement {
        start_end = exon_string.split("..").collect::<Vec<&str>>();
    }
    else{start_end = exon_string.split(|c| c=='(' ||c==')').collect::<Vec<&str>>()[1].split("..").collect::<Vec<&str>>(); }

    let start:i64 = match start_end[0].parse::<i64>(){
        Ok(val) => val,
        Err(e) => {eprintln!("Unparseable value of {} found in exon, error was {}", start_end[1], e);
        0},
    };
    let end:i64 = match start_end[1].parse::<i64>(){
        Ok(val) => val,
        Err(e) => {eprintln!("Unparseable value of {} found in exon, error was {}", start_end[1], e);
        0},
    };

    Exon{start, end, complement}
}

pub fn parse_taxonomy(first_line: &str, lines: &mut Peekable<Lines<BufReader<File>>>) -> Vec<String>{
    let mut taxonomy:Vec<String> = Vec::new();
    let first_fields = first_line.split(';').collect::<Vec<&str>>();
    assert!(first_line.starts_with("OC"), "Expected OC line for taxonomy, found {}", first_line);
    let special_case = first_fields[0].trim_start_matches("OC").trim();
    if special_case.split_whitespace().collect::<Vec<&str>>().len() > 1{
        taxonomy.push(special_case.replace(" ", "_").to_string());
    }
    else{
       taxonomy.push(special_case.to_string()); 
    }

    for field in first_fields[1..].iter(){
        let new_field= field.trim().replace(" ", "_");
         taxonomy.push(new_field.trim_end_matches('.').to_string());
    }
    let mut next_line = match lines.peek(){
        Some(Ok(line)) => line.to_string(),
        _ => panic!("Unexpected end of file")
    };
    while next_line.starts_with("OC"){
        let line = (&mut *lines).next().unwrap().unwrap();
        let fields = line.split(';').collect::<Vec<&str>>();
        let special_case = fields[0].trim_start_matches("OC").trim();
        if special_case.split_whitespace().collect::<Vec<&str>>().len() > 1{
            taxonomy.push(special_case.replace(" ", "_").to_string());
        }
        else{
            taxonomy.push(special_case.to_string()); 
        }
        for field in fields[1..].iter(){
            let new_field= field.trim().replace(" ", "_");
            taxonomy.push(new_field.trim_end_matches('.').to_string());
        }
        next_line = match lines.peek(){
            Some(Ok(line)) => line.to_string(),
            _ => panic!("Unexpected end of file")
        };
    }
    taxonomy
}
