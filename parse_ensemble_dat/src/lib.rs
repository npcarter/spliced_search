// parse_ensemble_dat: takes an ensemble dat file that describes a genome and
// parses it into DNA chunks that contain proteins along with vectors that
// describe their introns, exons, start points, and end points.
use std::fs::File;
use std::iter::Peekable;
use std::io::{BufReader, Lines};
use std::cmp::{min, max};
// For exons, need to know where they start and end in the contig,
// plus whether they are complementary to the DNA sequence.
pub struct Exon{
    pub start: i64,
    pub end: i64,
    pub complement: bool,
}

// For introns, need to know where they start and end in the contig.
// and whether they are complementary to the DNA sequence.
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
    // DNA sequence start
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

// translate a DNA sequence into a protein sequence
pub fn translate_sequence(mrna_sequence: &Vec<char>, codon_start: i8) -> Vec<char>{
    let mut protein_sequence:Vec<char> = Vec::new();
    let dna_sequence = mrna_sequence[codon_start as usize ..].to_vec();
   // if dna_sequence.len() % 3 != 0 { panic!("DNA sequence must be a multiple of three characters long"); }
    for i in 0..(dna_sequence.len()/3){
        protein_sequence.push(translate_codon(&dna_sequence[i*3..(i*3)+3].to_vec()));
    }

    protein_sequence
}




// Because the dat file gives the DNA sequence for the contig at the end of its
// description, after we've parsed the contig, we need to go through and fill in
// the protein sequences for each gene.  This function handles that.
pub fn parse_contig(contig: &mut Contig){
  //  println!("Contig {} has length {} and organism {}", contig.name, contig.length, contig.organism);
    for gene in contig.genes.iter_mut() {
        for protein in gene.proteins.iter_mut() {
            let mut check_sequence: Vec<char> = Vec::new();
            if !protein.exons[0].complement {
                for index in protein.mrna_start..protein.mrna_end + 1
                { protein.dna_sequence.push(contig.dna_sequence[index as usize - 1]); }

                /*     //need to handle first exon separately because coding sequence may not start at start of exon
                     for index in (protein.coding_start + protein.exons[0].start)..protein.exons[0].end+1{
                         check_sequence.push(protein.dna_sequence[(index -protein.mrna_start) as usize ]);
                     }
     */
                for exon in protein.exons.iter_mut() {
                    for index in max(protein.coding_start + protein.mrna_start, exon.start)..(min(exon.end, protein.coding_end + protein.mrna_start) + 1) {
                        check_sequence.push(protein.dna_sequence[(index - protein.mrna_start) as usize]);
                    }
                    //println!("{}, {}, length={}", exon.start-protein.mrna_start, min(exon.end, protein.coding_end+protein.mrna_start)-protein.mrna_start, check_sequence.len());
                }
            } else {
                for index in (protein.mrna_start..protein.mrna_end + 1).rev()
                {
                    protein.dna_sequence.push(complement(contig.dna_sequence[index as usize - 1]));
                }
                //println!("{}: ", protein.name);
                for exon in protein.exons.iter_mut() {
                    //println!("start: {}, end: {}", max(protein.coding_start + protein.mrna_start,exon.start), min(exon.end, protein.coding_end+ protein.mrna_start));
                    for index in (max(protein.coding_start + protein.mrna_start, exon.start)..(min(exon.end, protein.coding_end + protein.mrna_start) + 1)).rev() {
                        check_sequence.push(protein.dna_sequence[(protein.mrna_end- index)as usize]);
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
                for i in 1..translated_sequence.len() {
                    if translated_sequence[i] != protein.protein_sequence[i+start_offset] && translated_sequence[i] != 'X' && protein.protein_sequence[i+start_offset] != 'X' && !(translated_sequence[i] == '*' && protein.protein_sequence[i+start_offset] == 'U'){
                        found_missmatch = true;
                        println!("missmatch at position {}: {} vs {}", i, translated_sequence[i], protein.protein_sequence[i]);
                    }
                }
                if found_missmatch{
                    println!("Error: protein sequence {} does not match translated sequence {}", protein.protein_sequence.iter().collect::<String>(), translated_sequence.iter().collect::<String>());
                    let first_mismatch = translated_sequence[1..].iter()
                        .zip(protein.protein_sequence[1..].iter())
                        .position(|(a, b)| a != b);
                    println!("First mismatch at position {}.", first_mismatch.unwrap())
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
        introns: Vec::new(), exons: Vec::new(), codon_start: 0, dna_sequence: Vec::new(), protein_sequence: Vec::new()};

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

        for exon_string in line.split_whitespace().collect::<Vec<&str>>()[1].trim_end_matches(")").trim_end_matches(',').split(',') {
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

    let mut cds_exons = cds_fields[2].trim_start_matches("join(").trim_end_matches(")").trim_end_matches(',').split(',').collect::<Vec<&str>>();
    let first_cds_exon = parse_exon(cds_exons[0].trim_start_matches('(').trim_end_matches(')'));
    let mut last_cds_exon = parse_exon(cds_exons[cds_exons.len()-1].trim_start_matches('(').trim_end_matches(')'));

    // now, see if there are any other lines of exons
    cds_line = (&mut *lines).next().unwrap().unwrap();
    cds_fields = cds_line.split_whitespace().collect::<Vec<&str>>();
    while !cds_fields[1].starts_with("/"){
        cds_exons = cds_fields[1].trim_start_matches("join(").trim_end_matches(")").trim_end_matches(',').split(',').collect::<Vec<&str>>();
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


// generate an intron from the string that describes it.
pub fn parse_exon(exon_string: &str) -> Exon{
    let start_end:Vec<&str>;
    let complement:bool = exon_string.starts_with("complement");
    if !complement {
        start_end = exon_string.split("..").collect::<Vec<&str>>();
    }
    else{start_end = exon_string.split(|c| c=='(' ||c==')').collect::<Vec<&str>>()[1].split("..").collect::<Vec<&str>>(); }

    let start:i64 = start_end[0].parse::<i64>().unwrap();
    let end:i64 = start_end[1].parse::<i64>().unwrap();

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
