use std::fs::File;
use std::io::{BufReader};
use clap::Parser;
use noodles_gtf as gtf;
use noodles_fasta as fasta;
use noodles_gff as gff;
use noodles_core;
use parse_ncbi_gtf::*;
use bstr::{ByteSlice};
use std::collections::HashMap;
use std::str;
use rand::prelude::*;

// Configure command-line arguments
#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {
    basename: String,
}


fn main() {
    let args = Cli::parse();
    let protein_file:String;
    let mrna_file:String;
    let genome_file:String;
    let gtf_file: String;
    let mut count: u64 = 0;
    if args.basename.ends_with('/'){
        protein_file = args.basename.clone() + "protein.faa"; 
        mrna_file = args.basename.clone() + "rna.fna";
        genome_file = args.basename.clone() + "genomic.fna";
        gtf_file = args.basename.clone() + "genomic.gtf";
    }
    else{
        protein_file = args.basename.clone() + "/protein.faa"; 
        mrna_file = args.basename.clone() + "/rna.fna";
        genome_file = args.basename.clone() + "/genomic.fna";
        gtf_file = args.basename.clone() + "/genomic.gtf";
    }    


    // Three HashMaps we'll use to hold the protein, mrna, and raw sequence data from the various files NCBI provides
    let mut ncbi_mrnas:  HashMap<String, NcbiSequence> = HashMap::new();
    let mut ncbi_proteins: HashMap<String, NcbiSequence> = HashMap::new();
    let mut ncbi_genomes: HashMap<String, noodles_fasta::Record> = HashMap::new();
    // because the mrna and protein entries aren't in the same order, we have to go through 
    // both files and create a hashmap
    
 
     // First, read through the mrna file and create ncbi_sequence entries for all the mRNAs in it
    let file = File::open(mrna_file.clone());
    let mut mrna_reader : fasta::io::Reader<BufReader<File>>;
    match file {
        Ok(file) => {
            mrna_reader = fasta::io::Reader::new(BufReader::new(file));
        }
    
        Err(e) => {
            eprintln!("Error opening mrna file {} : {}", mrna_file, e);
            std::process::exit(1);
        }
    }
    for result in mrna_reader.records(){
        let mut the_entry = NcbiSequence::default();
        let record = result.expect("Error reading record from mrna file");

        let mrna_name: String = match str::from_utf8(record.name()){
            Ok(val) => {
                val.to_string()
            },
            Err (e) => {
                panic!("Couldn't convert mrna name to String, error was {}", e);
            }
        };

        the_entry.id = mrna_name.to_string();
        the_entry.mrna_sequence = str::from_utf8(record.sequence().as_ref()).expect("Invalid UTF-8 sequence in protein").chars().collect();
        if ncbi_mrnas.insert(mrna_name.to_string(), the_entry) != None{
            eprintln!("Duplicate mRNA key {} found", mrna_name);
        }
 
    }
    println!("{} coding mRNAs found", ncbi_mrnas.len());
    // now, add the proteins to the HashMap
    let file = File::open(protein_file.clone());
    let mut protein_reader : fasta::io::Reader<BufReader<File>>;
    match file {
        Ok(file) => {
            protein_reader = fasta::io::Reader::new(BufReader::new(file));
        }
    
        Err(e) => {
            eprintln!("Error opening protein file {} : {}", protein_file, e);
            std::process::exit(1);
        }
    }

    for result in protein_reader.records(){
        let record = result.expect("Error reading record from protein file");
        let mut the_entry = NcbiSequence::default();
        let id_string = match str::from_utf8(record.name()){
            Ok(val) => {
                val
            },
            Err (e) => {
                panic!("Couldn't convert protein name to String, error was {}", e);
            }
        };

        let sequence: Vec<char> = str::from_utf8(record.sequence().as_ref()).expect("Invalid UTF-8 sequence in protein").chars().collect();
        if sequence[0] != 'M' && id_string[0..2] == "NP_".to_string() {
            eprintln!("Non-methionine start AA of {} found in protein {}", sequence[0], id_string);
        }
        the_entry.id = id_string.to_string();
        the_entry.protein_sequence = sequence;
        if ncbi_proteins.insert(id_string.to_string(), the_entry) != None{
            eprintln!("Duplicate protein key {} found", id_string);
        }
    }   
    println!("{} proteins found", ncbi_proteins.len());

    // and build a name, Record HashMap for the genome pieces.
    // Keep these as noodles_fasta records because that package
    // has good functions for subsetting records.
    // This will bloat our memory usage, but avoids being dependent
    // on the GTF file referencing chunks of the genome in the order they appear in the file
    let file = File::open(genome_file.clone());
    let mut genome_reader : fasta::io::Reader<BufReader<File>>;
    match file {
        Ok(file) => {
            genome_reader = fasta::io::Reader::new(BufReader::new(file));
        }
    
        Err(e) => {
            eprintln!("Error opening genome file {} : {}", genome_file, e);
            std::process::exit(1);
        }
    }
    for result in genome_reader.records(){
        let record = result.expect("Error reading record from genome file");

        let id: String = match str::from_utf8(record.name()){
            Ok(val) => {
                val.to_string()
            },
            Err (e) => {
                panic!("Couldn't convert genome name to String, error was {}", e);
            }
        };
        if ncbi_genomes.contains_key(&id){
            eprintln!("Duplicate contig {} found", id);
        }
        ncbi_genomes.insert(id, record);
    }

    // Ok, we're finally ready to start parsing the GTF file.
    // Open the GTF file and create a reader
    let file = File::open(gtf_file.clone());
    let mut gtf_reader : gtf::io::Reader<BufReader<File>>;
    match file {
        Ok(file) => {
            gtf_reader = gtf::io::Reader::new(BufReader::new(file));
        }
        Err(e) => {
            eprintln!("Error opening GTF file {} : {}", gtf_file, e);
            std::process::exit(1);
        }
    }

    // Iterate through GTF records until we find the start of a gene.
    // If we're in a protein_coding gene, build a vector of the records that describe the gene.

    // When we find the start of a gene, if we're currently in a protein-coding gene, process that
    // into a gene record and handle it.  Otherwise, check whether the new gene is protein-coding
    // and make sure we have the FASTA entry that corresponds to it.

    let mut in_protein_coding = false;
    let mut gene_gtf_records = Vec::new();

    // If we get here, we've opened both the GTF and FASTA files successfully
    for record_result in gtf_reader.record_bufs() {             
        let the_recbuf = record_result.expect("Unable to get GTF record buffer from file");

        if the_recbuf.ty() != "gene" {
            if in_protein_coding{
                gene_gtf_records.push(the_recbuf.clone());
                }
            continue; // Skip non-gene records
            }
        else{  // This is a gene record, make sure it's for a protein coding gene
            if in_protein_coding{ // We just finished reading the records for a protein coding gene
                        // so process it
                let the_gene = parse_gtf_gene(&gene_gtf_records, &ncbi_genomes, &mut ncbi_proteins, & mut ncbi_mrnas);
                match the_gene {
                    Some(gene) => { count += write_gene_traindata(gene)},// Found a valid gene, write the training data for it.
                    None => {},
                }
              //  write_gene_traindata(the_gene);
            }
            gene_gtf_records = Vec::new(); // clear this, since we're starting a new gene.
            let attrib = the_recbuf.attributes();
            match attrib.get(b"gene_biotype"){
                None =>{
                    in_protein_coding = false;
                    continue; 
                }
                Some(val) =>{
                    let the_value = val;
                    if the_value.as_string().expect("Couldn't get gene_biotype as string") != "protein_coding"{
                        in_protein_coding = false;
                        continue;
                    }
                    else{
                        in_protein_coding = true;
                        gene_gtf_records.push(the_recbuf.clone());
                    }
                }
            }
                
            // If we make it here without hitting a continue, the GTF record is the start of a gene with the protein_coding type,
            // so start building the vector of records that represent the gene


                    
        }   
    }
    println!("Found {} usable proteins", count);
}   



// parse the GTF records that describe a gene into a structure
// requires that the gene is a protein coding gene
fn parse_gtf_gene(gtf_records:&Vec<gff::feature::RecordBuf> , ncbi_genomes:&HashMap<String, noodles_fasta::Record>, ncbi_proteins:& mut HashMap<String, NcbiSequence>, 
    ncbi_mrnas:&HashMap<String, NcbiSequence>) -> Option<Gene>{
    let mut the_gene = Gene{
        gene_id: String::new(),
        start: 0,
        end: 0,
        strand: Strand::Unknown,
        proteins: Vec::new(),
    };

    let mut the_protein = Protein::default();
    let mut rng = rand::rng(); 
    for record in gtf_records {
        match record.ty().to_str() { // what type of record is this
            Ok("gene")  => { // start of new gene, create the record
                assert!(the_gene.gene_id == ""); // There should only be one gene record in a gene.

                let attrib = record.attributes();
                match attrib.get(b"gene_id") {
                    None => {
                        eprintln!("GTF gene entry gene_id attribute wasn't found");
                        the_gene.gene_id = "NONE".to_string();
                    }
                    Some(val) => {
                        assert!(attrib.get(b"gene_biotype").expect("Parse_gtf_gene couldn't get gene_biotype of gene record").as_string().expect("Couldn't get gene_biotype as string") == "protein_coding");
                        the_gene.gene_id = val.as_string().expect("GTF gene entry gene_id attribute wasn't a string").to_string();
                    }
                }
                match record.strand(){ // This is a bit ugly because of duelling Strand types between this code and noodles_gff
                    noodles_gff::feature::record::Strand::None => {                 
                        eprintln!("Gene record found with no strand");
                        std::process::exit(1);
                    }
                    noodles_gff::feature::record::Strand::Unknown => {                 
                        eprintln!("Gene record found with unknown strand");
                        std::process::exit(1);
                    }
                    noodles_gff::feature::record::Strand::Forward =>{
                        assert!(the_gene.strand!= Strand::Reverse); //Complain if gene seems to switch strands
                        the_gene.strand = Strand::Forward;
                    }
                    noodles_gff::feature::record::Strand::Reverse =>{
                        assert!(the_gene.strand != Strand::Forward); // Complain if gene seems to switch strands
                        the_gene.strand = Strand::Reverse;
                    }
                } 
            }
            
            Ok("transcript") => {
            // start of new protein, so push old one onto gene list if there is one.
                if the_protein.name != None && sanity_check(&mut  the_protein){ // There was a previous protein, and it was well-defined
                    match ncbi_proteins.get_mut(&the_protein.name.clone().unwrap()){
                        None => {
                            eprintln!("NCBI protein entry for {} not found, skipping", the_protein.name.clone().unwrap());
                        }
                        Some(val) => {
                            if val.already_found{
                        //        eprintln!("NCBI protein entry for {} already found, skipping", the_protein.name.clone().unwrap());
                            }
                            else{
                                val.already_found = true;
                                the_gene.proteins.push(the_protein);
                            }
                        }
                    }
                };
                the_protein = Protein::default();
                // Grab what information we can from this record
                let fasta_source: String = record.reference_sequence_name().to_str().expect("Couldn't convert GTF source to string").to_string();
                let fasta_record = ncbi_genomes.get(&fasta_source).expect("Couldn't find DNA for source"); 
                the_protein.mrna_start = record.start().get() ;
                the_protein.mrna_end = record.end().get();

                let seq_length = the_protein.mrna_end - the_protein.mrna_start +1;
                let start_padding = std::cmp::max(rng.random_range(0..seq_length/10), 100);
                let end_padding = std::cmp::max(rng.random_range(0..seq_length/10), 100);

                the_protein.mrna_start = std::cmp::max(the_protein.mrna_start.saturating_sub(start_padding), 1);
                the_protein.mrna_end = std::cmp::min(the_protein.mrna_end.saturating_add(end_padding), fasta_record.sequence().len());
                let start = noodles_core::Position::try_from(the_protein.mrna_start).expect("Couldn't generate Position from mrna_start");
                let end = noodles_core::Position::try_from(the_protein.mrna_end).expect("Couldn't generate position from mrna_end");


                // Fetch the DNA sequence from the FASTA record
                match the_gene.strand{
                    Strand::Unknown => {                 
                        eprintln!("Gene record had unknown strand after transcript line.");
                        std::process::exit(1);
                    }
                    Strand::Forward => {
                        for residue in Vec::from(fasta_record.sequence().get(start..=end).expect("Couldn't extract DNA sequence from forward strand transcript")){
                            for c in (residue as char).to_uppercase(){ // this bit of ugliness brought to you by the fact 
                                // that to_uppercase() sometimes returns multiple characters, though that should never happen here.
                                the_protein.dna_sequence.push(c);
                            }  
                        }
                    }
                    Strand::Reverse => {
                        the_protein.dna_sequence = Vec::new(); // make sure this starts empty
                        // complement function appears to complement residues, but not reverse order.
                        for residue in fasta_record.sequence().slice(start..=end).expect("Couldn't extract DNA sequence from reverse strand transcript").complement().rev(){
                            match residue{
                                Err(e) => {
                                    eprintln!("Error encountered complementing DNA sequence: {}", e);
                                    std::process::exit(1);
                                }
                                Ok(val) => { for c in (val as char).to_uppercase(){ // this bit of ugliness brought to you by the fact 
                                    // that to_uppercase() sometimes returns multiple characters, though that should never happen here.
                                    the_protein.dna_sequence.push(c);
                                }}
                            }
                        }
                    }
                }
            }
            
            Ok("exon") => {  // parse the record and add an exon to the protein's list

                assert!(record.start().get() <= record.end().get()); // check to make sure our expectations about start and end are preserved
                // convert from 1-indexed positions within the FASTA record 
                // to 0-indexed offsets within the transcript
                // if reverse strand, handle the fact that we've reversed the transcript as part of complementing it
                if the_protein.exons.len() == 0{ // This is the first exon, get the NCBI mRNA,
                    let attrib = record.attributes();
                    let mut transcript:String = match attrib.get(b"transcript_id") {
                        None => {
                            eprintln!("Transcript ID attribute wasn't found in CDS record");
                            "bad_bad".to_string()
                        }
                        Some(val) => {
                            val.as_string().expect("GTF protein_id attribute wasn't a string").to_string()
                        }
                    };
                    if transcript.matches('_').count() > 1{ // variant naming
         //               print!("Variant transcript name {} found,", transcript);
                        let fields: Vec<&str> = transcript.split('_').collect();
                        transcript = fields[..fields.len()-1].join("_");
                    }
                    let the_entry = ncbi_mrnas.get(&transcript);
                    match the_entry{
                        None => {
                        }
                        Some(val) => {the_protein.ncbi_mrna = val.mrna_sequence.clone()}
                    }                      
                }
                let start:usize = match the_gene.strand{  
                    Strand::Forward => {
                        record.start().get() - the_protein.mrna_start
                    }
                    Strand::Reverse => {
                        the_protein.mrna_end -record.end().get()
                    }
                    Strand::Unknown => {
                        panic!("Reached impossible branch in computing exon start");
                    }
                };
                let end:usize = match the_gene.strand{
                    Strand::Forward => {
                        record.end().get() - the_protein.mrna_start
                    }
                    Strand::Reverse => {
                        the_protein.mrna_end -record.start().get()
                    }
                    Strand::Unknown => {
                        panic!("Reached impossible branch in computing exon start");
                    }
                };
                    record.start().get();

                let the_exon = Exon{
                    start: start,
                    end: end,
                    start_offset: None,
                };
                the_protein.exons.push(the_exon);
            }
            
            Ok("CDS") => { // coding segment.  
                // If this is the first coding segment, we need to extract the
                // name of the protein and the start position of the start codon
                if the_protein.name == None{  // This is the first coding segment of the protein
                    
                    // compute the offset from the start of the first codon to the first
                    // nucleotide in the start codon, accounting for reverse-strand effects
                    let start_codon_start = match the_gene.strand{
                        Strand::Forward => {
                            assert!(record.start().get() >= the_protein.mrna_start);  // codon start shouldn't be before mrna start
                            record.start().get() - the_protein.mrna_start
                        }
                        Strand::Reverse => {
                            assert!(record.end().get() <= the_protein.mrna_end); // codon start shouldn't be after mrna end 
                            // on reverse strand
                            the_protein.mrna_end - record.end().get()
                        }
                        Strand::Unknown => {
                            panic!("Reached impossible branch in computing start codon offset");
                        }
                    };
                    the_protein.coding_start = start_codon_start;
                    the_protein.start_codon = vec!(start_codon_start, start_codon_start+1, start_codon_start+2);
                    // retrieve the protein name from the CDS record
                    let attrib = record.attributes();
                    match attrib.get(b"protein_id") {
                        None => {
                            eprintln!("Protein_id attribute wasn't found in CDS record");
                            the_protein.name = Some("NONE".to_string());
                        }
                        Some(val) => {
                            the_protein.name = Some(val.as_string().expect("GTF protein_id attribute wasn't a string").to_string());
                        }
                    }
                    
                    let ncbi_protein_option = ncbi_proteins.get(&the_protein.name.clone().unwrap());
                    let ncbi_protein: & NcbiSequence= match ncbi_protein_option{
                        Some(val) => val,
                        None => {
                            eprintln!("Unable to find NCBI protein entry for {}", the_protein.name.clone().unwrap());
                            & NcbiSequence::default()
                        }
                    };
                    the_protein.ncbi_protein = ncbi_protein.protein_sequence.clone();
                    let attrib = record.attributes();
                    match attrib.get(b"exception"){
                        None =>{},                
                        Some(val) =>{
                            let the_value = val;
                            if the_value.as_string().expect("Couldn't get exception value as string") == b"alternative start codon"{
                                the_protein.alternative_start_codon = true;  // NCBI knows this has an alt start codon
                            }
                        }
                    }
                    match attrib.get(b"note"){
                        None => {},
                        Some (val) => {
                            let the_value = val;
                            let val_string = the_value.as_string().expect("Couldn't get note value as string").to_string();
                            if val_string.contains("non-AUG"){
                                the_protein.alternative_start_codon = true;  // NCBI knows this has an alt start codon
                            }
                            if val_string.contains("substitution")|| val_string.contains("aligns at")|| val_string.contains("frameshift"){
                                the_protein.wrong_transcript = true;  // NCBI knows this has an alt start codon
                            }
                            if val_string.contains("recoded as selenocysteine"){
                                the_protein.selenocystine = true;
                            }
                        }
                    }
                    match attrib.get(b"partial"){
                        None =>{},                
                        Some(val) =>{
                            let the_value = val;
                            if the_value.as_string().expect("Couldn't get exception value as string") == b"true"{
                                the_protein.wrong_transcript = true;  // NCBI knows this has an alt start codon
                            }
                            else{
                                eprintln!("Found partial with value {} in protein {:?}", the_value.as_string().expect("Couldn't get exception value as string"), the_protein.name);
                            }
                        }
                    }
                    match attrib.get(b"transl_except"){
                        None =>{},                
                        Some(val) =>{
                            let the_values = val;
                            for value in the_values{
                            let val_string = value.to_str().expect("Couldn't convert transl_except value to string").to_string();
                            if val_string.contains("aa:Sec"){ // other way of indicating selenocystine
                                the_protein.selenocystine = true;  // NCBI knows this has an alt start codon
                            }
                        }
                        }
                    }
                }
                match record.phase(){
                    Some(val) => {
                        assert!(the_protein.exons.len() >0); //We should have at least one exon by the time we hit a CDS entry 
                        let last_exon_index = the_protein.exons.len() -1;
                        match val {
                            noodles_gff::feature::record::Phase::Zero => {the_protein.exons[last_exon_index].start_offset = Some(CodonPosition::First);}
                            noodles_gff::feature::record::Phase::One => {the_protein.exons[last_exon_index].start_offset = Some(CodonPosition::Second);} 
                            noodles_gff::feature::record::Phase::Two => {the_protein.exons[last_exon_index].start_offset = Some(CodonPosition::Third);}
                        }
                    }
                    None => {
                        panic!("CDS entry found without phase");
                    }
                }
                
                // update the protein's coding_end value if necessary 
                let coding_end = match the_gene.strand{
                    Strand::Forward => {
                        assert!(record.end().get() >= the_protein.mrna_start);  // codon start shouldn't be before mrna start
                        std::cmp::max(record.end().get() - the_protein.mrna_start, the_protein.coding_end)
                    }
                    Strand::Reverse => {
                        assert!(record.start().get() <= the_protein.mrna_end); // codon start shouldn't be after mrna end 
                        // on reverse strand
                        std::cmp::max(the_protein.mrna_end - record.start().get(), the_protein.coding_end)
                    }
                    Strand::Unknown => {
                        panic!("Reached impossible branch in computing stop codon offset");
                    }
                };
                the_protein.coding_end = coding_end;
                // fill in default stop codon positions
                assert!(the_protein.coding_end + 3 < the_protein.dna_sequence.len()); // make sure we won't go off the end of the protein
                the_protein.stop_codon[0] = coding_end+1;  // This doesn't have to be handled differently for reverse strand because we've already reversed
                the_protein.stop_codon[1] = coding_end+2;
                the_protein.stop_codon[2] = coding_end+3;
            }

            Ok("start_codon") => {  // If the protein has one of these entries, check that 
                // it matches what we computed
                let start_codon_start = match the_gene.strand{
                    Strand::Forward => {
                        assert!(record.start().get() >= the_protein.mrna_start);  // codon start shouldn't be before mrna start
                        record.start().get() - the_protein.mrna_start
                    }
                    Strand::Reverse => {
                        assert!(record.end().get() <= the_protein.mrna_end); // codon start shouldn't be after mrna end 
                        // on reverse strand
                        the_protein.mrna_end - record.end().get()
                    }
                    Strand::Unknown => {
                        panic!("Reached impossible branch in computing start codon offset");
                    }
                };

                if record.end().get() - record.start().get() == 2{ // this is a well-formed start_codon entry, so check it
                    assert!(start_codon_start == the_protein.coding_start);
                }

                assert!((record.end().get() - record.start().get())+ 1 + the_protein.start_positions_found <=3); // make sure we haven't found too many start codon positions
                match the_gene.strand{
                    Strand::Forward => {
                        for i in record.start().get()..=record.end().get(){
                            the_protein.start_codon[the_protein.start_positions_found] = i - the_protein.mrna_start;
                            the_protein.start_positions_found+=1; 
                        }
                    },
                    Strand::Reverse => {
                        for i in (record.start().get()..=record.end().get()).rev(){
                            the_protein.start_codon[the_protein.start_positions_found] = the_protein.mrna_end - i;
                            the_protein.start_positions_found +=1;
                        }
                    },
                    Strand::Unknown => {
                        panic!("Reached impossible branch in computing start codon location");
                    }
                }
            //    else{
              //      println!("miss-formed start_codon entry found for protein {:?}", the_protein.name);
               // }
            }
            Ok("stop_codon") => {  // If the protein has one of these entries, check that 
                // it matches what we computed
                let stop_codon_start = match the_gene.strand{
                    Strand::Forward => {
                        assert!(record.end().get() <= the_protein.mrna_end);  // stop codon shouldn't be outside mrna
                        record.start().get() - the_protein.mrna_start
                    }
                    Strand::Reverse => {
                        assert!(record.start().get() >= the_protein.mrna_start); // codon start shouldn't be after mrna end 
                        // on reverse strand
                        the_protein.mrna_end - record.end().get()
                    }
                    Strand::Unknown => {
                        panic!("Reached impossible branch in computing start codon offset");
                    }
                };

                assert!((record.end().get() - record.start().get())+ 1 + the_protein.stop_positions_found <=3); // make sure we haven't found too many stop codon positions
                match the_gene.strand{
                    Strand::Forward => {
                        for i in record.start().get()..=record.end().get(){
                            the_protein.stop_codon[the_protein.stop_positions_found] = i - the_protein.mrna_start;
                            the_protein.stop_positions_found+=1; 
                        }
                    },
                    Strand::Reverse => {
                        for i in (record.start().get()..=record.end().get()).rev(){
                            the_protein.stop_codon[the_protein.stop_positions_found] = the_protein.mrna_end - i;
                            the_protein.stop_positions_found +=1;
                        }
                    },
                    Strand::Unknown => {
                        panic!("Reached impossible branch in computing start codon location");
                    }
                }
            //    else{
              //      println!("miss-formed start_codon entry found for protein {:?}", the_protein.name);
               // }
            }

            Ok("Selenocysteine") => {
                // special entry that records the presence of a selenocystine.  We don't care about that, so do nothing.
            }
            Err(e)=> {
                eprintln!("Unable to parse record type as string: {}", e);
                std::process::exit(1);
            }

            Ok(val) => {
                eprintln!("Unexpected record type field found: {}", val);
                std::process::exit(1);
            } 
        }
    }

    if sanity_check(&mut the_protein){
        match ncbi_proteins.get_mut(&the_protein.name.clone().unwrap()){
            None => {
                eprintln!("NCBI protein entry for {} not found, skipping", the_protein.name.clone().unwrap());
            }
            Some(val) => {
                if val.already_found{
                   // eprintln!("2 NCBI protein entry for {} already found, skipping", the_protein.name.clone().unwrap());
                }
                else{
                    val.already_found = true;
                    the_gene.proteins.push(the_protein);
                }
            }
        }
    }
    if the_gene.proteins.len() >0 {
        Some(the_gene)
    }
    else{
        None
    }
}

fn write_gene_traindata(the_gene: Gene)-> u64{
    let mut count = 0;
    for protein in the_gene.proteins{
 

        count +=1;

     //   println!("{}{}{}", protein.dna_sequence[protein.coding_end -2], protein.dna_sequence[protein.coding_end-1], protein.dna_sequence[protein.coding_end]);
    }
    count
}

// check a protein data structure to make sure it's been well-defined
// NCBI GTF seems to contain some broken proteins, so need to check
fn sanity_check(the_protein:&mut Protein)-> bool{
    match &the_protein.name{ // Only want curated proteins in training data
        None => {}, // This is a failure case, but don't fail here because we want to figure out what went wrong.
        Some(name) => {
            if name[0..3] == "YP_".to_string(){ // this protein had no associated transcript, ignore it
               return false;
            }
        }
    }

    if the_protein.dna_sequence.len() == 0 || the_protein.ncbi_mrna.len() == 0 || the_protein.ncbi_protein.len() == 0{
        // one or more of these weren't found 
        eprintln!("Protein {:?} had at least one bad sequence field. dna_sequence length was {}, ncbi_mrna length was {}, ncbi_protein length was {}", the_protein.name, the_protein.dna_sequence.len(), 
    the_protein.ncbi_mrna.len(), the_protein.ncbi_protein.len());
        return false;
    }
    
    if the_protein.coding_start == std::usize::MAX || the_protein.coding_end == 0 || the_protein.coding_end <= the_protein.coding_start {
        eprintln!("Protein {:?} had invalid coding start and/or end", the_protein.name);
        // These fields weren't set or were set wrong
        return false;
    }
    if the_protein.mrna_start == std::usize::MAX || the_protein.mrna_end == 0 || the_protein.mrna_end <= the_protein.mrna_start {
        eprintln!("Protein {:?} had invalid coding start and/or end", the_protein.name);
        // These fields weren't set or were set wrong
        return false;
    }

    if the_protein.exons.len() == 0{  // poorly-defined protein had no exons
        eprintln!("protein had no exons");
        return false;
    }
    
    if the_protein.name == None{ // This field wasn't set
        eprintln!("Protein had missing name, rejecting");
        return false;
    }
    if the_protein.wrong_transcript{
      //  eprintln!("Rejecting {:?} because GTF noted that it did not match the reference transcript", the_protein.name);
        return false;
    }
    if the_protein.alternative_start_codon  == false && (the_protein.dna_sequence[the_protein.start_codon[1]] != 'T' || 
        the_protein.dna_sequence[the_protein.start_codon[2]] != 'G')
    {
        println!("Rejecting {:?} for bad start codon {}{}{}", the_protein.name, the_protein.dna_sequence[the_protein.start_codon[0]] as char, 
            the_protein.dna_sequence[the_protein.start_codon[1]] as char,
            the_protein.dna_sequence[the_protein.start_codon[2]] as char);
        return false;
    }
    // if we get this far, protein has passed the basic checks, see if our computed mRNA matches the NCBI one and then try to hand translate
    for exon in &the_protein.exons{
        for i in std::cmp::max(exon.start, the_protein.coding_start)..=std::cmp::min(exon.end, the_protein.coding_end){
            the_protein.computed_coding.push(the_protein.dna_sequence[i]);
        }
    }
    if the_protein.computed_coding.len() %3 != 0{
        eprintln!("Rejecting {:?} because computed coding sequence length {} was not a multiple of three", the_protein.name, 
        the_protein.computed_coding.len());
        return false;
    }
    
    if the_protein.computed_coding.len() /3 != the_protein.ncbi_protein.len(){
        eprintln!("Rejecting {:?} because computed protein length {} did not match reference length {}", the_protein.name, 
        the_protein.computed_coding.len()/3, the_protein.ncbi_protein.len());
        return false;
    }
 
    // Don't check the translation of the first AA if the protein has an alternative start codon, as it won't match
    let start_check = match the_protein.alternative_start_codon{
        true => 1,
        false => 0,
    };

    for i in (0..the_protein.computed_coding.len()).step_by(3){
        the_protein.computed_protein.push(translate_codon(the_protein.computed_coding[i..i+3].to_vec(), the_protein.selenocystine));
    }

    for i in start_check..std::cmp::min(the_protein.computed_protein.len(), the_protein.ncbi_protein.len()){
        if the_protein.computed_protein[i] != the_protein.ncbi_protein[i]{
eprintln!("Rejecting {:?} because computed protein did not match reference protein at position {}, {} vs {}", the_protein.name, i, 
      the_protein.computed_protein[i], the_protein.ncbi_protein[i]);
            return false;
        }
    }

    // check the stop codon
    let stop_codon: String = vec!(the_protein.dna_sequence[the_protein.stop_codon[0]], 
        the_protein.dna_sequence[the_protein.stop_codon[1]], the_protein.dna_sequence[the_protein.stop_codon[2]]).into_iter().collect();

        let valid_stop: bool = match stop_codon.as_str(){
            "TAA" => {true},
            "TAG" => {true},
            "TGA" => {true},
            _ => {false},
        };
        if !valid_stop{
            eprintln!("Rejecting {:?} because of invalid stop codon {}", the_protein.name, stop_codon);
            return false;
        }
    // If we get this far, we've passed all the checks, so accept the protein
    true
}