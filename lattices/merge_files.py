from genus import write_header_to_file

def merge_files(fnames, out_fname):
    write_header_to_file(out_fname)
    for fname in fnames:
        with open(fname) as f:
        # only needed if they have headers
        #    lines = f.readlines()[3:]
            lines = f.readlines()
        with open(out_fname, "a") as f:
            f.writelines(lines)
