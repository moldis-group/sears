import sys
  
def replace_N_with_H_inplace(xyz_path):
    with open(xyz_path, "r") as f:
        lines = f.readlines()

    if len(lines) < 3:
        print("File too short to be valid XYZ.")
        return

    # First two lines: atom count + comment
    header = lines[:2]
    body = lines[2:]

    new_body = []
    for line in body:
        parts = line.strip().split()
        if len(parts) < 4:
            new_body.append(line)
            continue

        atom = parts[0]
        coords = parts[1:]

        if atom == "N":
            atom = "H"

        new_line = f"{atom:2s}  {coords[0]}  {coords[1]}  {coords[2]}\n"
        new_body.append(new_line)

    # Write back to same file (in place)
    with open(xyz_path, "w") as f:
        f.writelines(header)
        f.writelines(new_body)

    print(f"Replaced N ?~F~R H in: {xyz_path}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python replace_N_with_H_inplace.py file.xyz")
        sys.exit(1)

    xyz_file = sys.argv[1]
    replace_N_with_H_inplace(xyz_file)
