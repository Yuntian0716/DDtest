import argparse
import csv
import h5py
import os

def load_barcodes_from_h5(filepath):
    # Open the HDF5 file and extract the barcodes stored in the "matrix/barcodes" dataset.
    with h5py.File(filepath, "r") as f:
        barcodes = f["matrix"]["barcodes"][:]
        return [barcode.decode("utf-8") if isinstance(barcode, bytes) else barcode for barcode in barcodes]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--filtered-h5', required=True, help='Path to filtered_feature_bc_matrix.h5')
    parser.add_argument('--raw-h5', required=True, help='Path to raw_feature_bc_matrix.h5')
    parser.add_argument('--output', default="singlecell.csv", help='Output CSV file path or directory')
    args = parser.parse_args()

    # Load barcodes from the filtered and raw files.
    pass_qc_barcodes = set(load_barcodes_from_h5(args.filtered_h5))
    all_barcodes = load_barcodes_from_h5(args.raw_h5)

    # Sanity check: Ensure that all filtered barcodes are present in the raw barcode list.
    if not pass_qc_barcodes.issubset(set(all_barcodes)):
        raise ValueError("Some pass QC barcodes are not in the full barcode list?")

    # Check if output is a directory; if so, append the default filename.
    output_path = args.output
    if os.path.isdir(output_path):
        output_path = os.path.join(output_path, "singlecell.csv")

    # Write the output CSV with the header.
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["barcode", "is__cell_barcode"])
        for barcode in all_barcodes:
            writer.writerow([barcode, int(barcode in pass_qc_barcodes)])

if __name__ == "__main__":
    main()