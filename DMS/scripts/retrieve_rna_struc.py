import subprocess as sp
import glob
from pdf2image import convert_from_path
import os
from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw
import pickle
from pathlib import Path

def get_DMS_values(file_name, clamp_value = 0.2):
    
    file_ = open(file_name)
    file_lines = file_.readlines()
    value_list = file_lines[1].split(" ")[:-1]
    value_list = [str(min(float(v),clamp_value)) + "\n" for v in value_list]

    return value_list


def color_ps_file(file_name, dms_vaues_path, function_file):

    with open(function_file) as f:
        func_lines = f.readlines()

    try:
        replace_idx = next(i for i, ln in enumerate(func_lines) if ln == "REPLACE_ME_WITH_VALUES\n")
    except StopIteration:
        raise ValueError("Placeholder line 'REPLACE_ME_WITH_VALUES' not found in function_file.")

    values = get_DMS_values(dms_vaues_path)
    func_lines[replace_idx:replace_idx + 1] = values

    with open(file_name) as f:
        ps_lines = f.readlines()

    try:
        at_idx = next(i for i, ln in enumerate(ps_lines)
                      if ln == "% switch off outline pairs or bases by removing these lines\n")
    except StopIteration:
        raise ValueError("Line for insertion not found in PS file.")

    try:
        draw_idx = next(i for i, ln in enumerate(ps_lines) if ln == "drawoutline\n")
    except StopIteration:
        raise ValueError("'drawoutline' line not found in PS file.")

    ps_lines[at_idx:at_idx] = func_lines
    draw_insert_at = draw_idx + len(func_lines)

    draw_list = ["/invert false def\n", "drawreliability\n", "0.1 0.1 colorbar\n"]
    ps_lines[draw_insert_at:draw_insert_at] = draw_list

    with open(file_name, "w") as f:
        f.writelines(ps_lines)

    return


def RNAfold(prodigal_cmd, fasta_file, constraint):
    
    fasta_file_preffix = fasta_file.rsplit('.', 1)[0]
    output_pdf = fasta_file_preffix + '_proteins.fa'
    log_file = fasta_file_preffix + '_RNAfold.log'
    prodigal_cmd += ' {input_fasta} --filename-full'
    if constraint == True: prodigal_cmd += " -C" 
        
    prodigal_cmd = prodigal_cmd.format(prodigal=prodigal_cmd, input_fasta=fasta_file)

    with open(log_file, 'w') as lf:
        sp.call(prodigal_cmd.split(), stdout=lf)
        
    return


def get_concat_v(im1, im2, im3):
    dst = Image.new('RGB', (im1.width, im1.height + im2.height + im3.height))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (0, im1.height))
    dst.paste(im3, (0, im1.height + im2.height))
    return dst



def expand_white_space(im1):
    dst = Image.new('RGB', (im1.width, im1.height + 80), color = "white")
    dst.paste(im1, (0, 80))
    
    return dst
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser("retrieve RNA seq structure")
    parser.add_argument("-f", "--function_file", required=True, help="path to function file", type=str)
    parser.add_argument("-t", "--target_folder", required=True, help="output folder for PS/JPG files", type=str)
    parser.add_argument("-s", "--structure_folder", required=True, help="folder containing input FASTA/constraint files", type=str)
    parser.add_argument("-a", "--summary", required=True, help="summary output file path", type=str)
    parser.add_argument("-z", "--type", required=True, help="font file for labels (e.g., a .ttf path)", type=str)
    parser.add_argument("-d", "--dms_path", required=True, help="path to saved dms files", type=str)
    parser.add_argument("-u", "--untreated_path", required=True, help="path to saved untreated files", type=str)

    args, _ = parser.parse_known_args()

    in_path = Path(args.target_folder)
    in_path.mkdir(parents=True, exist_ok=True)

    fasta_files_to_run = sorted(glob.glob(str(Path(args.structure_folder) / "*")))

    dms_path = args.dms_path
    untreated_val_path = args.untreated_path



    for constraint_file in fasta_files_to_run:
        cpath = Path(constraint_file)
        basename = cpath.name
        sub_string = basename.split("_")[0]

        dms_file_treated = basename.split("treated")[0] + "treated.fa"
        dms_file_untreated = basename.split("treated")[0] + "untreated.fa"

        with open(cpath) as fh:
            lines = fh.readlines()
        id_only, sequence_only, constraint_only = lines[0], lines[1], lines[2]

        constraint_file_mod_path = Path("./") / (cpath.stem + "_mod.fa")
        with open(constraint_file_mod_path, "w") as fh:
            fh.write(id_only)
            fh.write(sequence_only)

        untreated_file = (
            basename.split("treated")[0]
            + "untreated_"
            + str(untreated_val_path).split("_")[-1].split("/")[0]
            + ".fa"
        )
        untreated_path = Path(untreated_val_path) / untreated_file

        RNAfold("RNAfold", str(constraint_file_mod_path), False)
        current_files = sorted(Path(".").glob("*.ps"))
        match = [p for p in current_files if sub_string in p.name]
        if not match:
            raise FileNotFoundError(f"No .ps produced for {constraint_file_mod_path} (sequence-only)")
        sequence_name = in_path / (match[0].name.split("ss")[0] + "_only_seq_ss.ps")
        os.rename(str(match[0]), str(sequence_name))

        RNAfold("RNAfold", str(cpath), True)
        current_files = sorted(Path(".").glob("*.ps"))
        match = [p for p in current_files if sub_string in p.name]
        if not match:
            raise FileNotFoundError(f"No .ps produced for {cpath} (treated)")
        treated_name = in_path / (match[0].name.split("ss")[0] + "_treated_ss.ps")
        os.rename(str(match[0]), str(treated_name))
        color_ps_file(
            str(treated_name),
            str(Path(dms_path) / (dms_file_treated.split("_constraints")[0] + "_treated.fa")),
            args.function_file,
        )

        RNAfold("RNAfold", str(untreated_path), True)
        current_files = sorted(Path(".").glob("*.ps"))
        match = [p for p in current_files if sub_string in p.name]
        if not match:
            raise FileNotFoundError(f"No .ps produced for {untreated_path} (untreated)")
        untreated_name = in_path / (match[0].name.split("ss")[0] + "_untreated_ss.ps")
        os.rename(str(match[0]), str(untreated_name))
        font = ImageFont.truetype(args.type, size=30)
        color_ps_file(
            str(untreated_name),
            str(Path(dms_path) / (dms_file_untreated.split("_constraints")[0] + "_untreated.fa")),
            args.function_file,
        )

        # cleanup temp file
        try:
            os.remove(constraint_file_mod_path)
        except FileNotFoundError:
            pass

        psimage1 = expand_white_space(Image.open(sequence_name))
        draw1 = ImageDraw.Draw(psimage1)
        title = basename.split("_constraints")[0]
        l1 = draw1.textlength(title, font=font)
        draw1.text((psimage1.width / 2 - l1 / 2, 5), title, (0, 0, 0), font=font)
        l1 = draw1.textlength("sequence only", font=font)
        draw1.text((psimage1.width / 2 - l1 / 2, 40), "sequence only", (0, 0, 0), font=font)

        psimage2 = expand_white_space(Image.open(treated_name))
        draw2 = ImageDraw.Draw(psimage2)
        l2 = draw2.textlength("1lbeta-treated", font=font)
        draw2.text((psimage2.width / 2 - l2 / 2, 30), "1lbeta-treated", (0, 0, 0), font=font)

        psimage3 = expand_white_space(Image.open(untreated_name))
        draw3 = ImageDraw.Draw(psimage3)
        l3 = draw3.textlength("untreated", font=font)
        draw3.text((psimage3.width / 2 - l3 / 2, 30), "untreated", (0, 0, 0), font=font)

        dst = get_concat_v(psimage1, psimage2, psimage3)
        out_jpg = in_path / (basename.split("treated")[0] + ".jpg")
        dst.save(out_jpg)

        for p in (sequence_name, treated_name, untreated_name):
            try:
                os.remove(p)
            except FileNotFoundError:
                pass

        for log in glob.glob("./*.log"):
            try:
                os.remove(log)
            except OSError:
                print("already gone")

    with open(args.summary, "w") as sf:
        sf.write("finished")

