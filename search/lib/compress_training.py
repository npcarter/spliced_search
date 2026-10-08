from dataloader import compress
import argparse

def main(input_file_path, output_file_path, sequence_chunk_size, name_length, validate, test):
    compress.compress_file(input_file_path, output_file_path, sequence_chunk_size, name_length, 1024, validate, test)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Compresses text-format training arguments into compressed/Zarr format")

    parser.add_argument("input_file_path")
    parser.add_argument("output_file_path")
    parser.add_argument("--chunk", required=False, default=16384, type=int)
    parser.add_argument("--nl", required=False, default=40, type=int)
    parser.add_argument("--validate", required=False, default=0.10, type=float)
    parser.add_argument("--test", required=False, default=0.01)


    args = parser.parse_args()
    
    main(args.input_file_path, args.output_file_path, args.chunk, args.nl, args.validate, args.test)