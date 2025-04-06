import subprocess as sp
import glob
from pdf2image import convert_from_path
import os
from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw

#CCT2-NM_006431

def get_DMS_values(file_name, clamp_value = 0.2):
    
    file_ = open(file_name)
    file_lines = file_.readlines()
    
    value_list = file_lines[1].split(" ")[:-1]
    value_list = [str(min(float(v),clamp_value)) + "\n" for v in value_list]

    return value_list


def color_ps_file(file_name, dms_vaues_path, function_file):

    file2_ = open(function_file)
    lines = file2_.readlines()
    file2_.close()
    replace_int = None
    

    for enum,line in enumerate(lines):

        if str(line) == "REPLACE_ME_WITH_VALUES\n":
            replace_int = enum
            
            
    values = get_DMS_values(dms_vaues_path)
    print(values)
    
    lines[replace_int:replace_int+1] = values
    
    
    file_ = open(file_name)
    lines_read = file_.readlines()
    at_line = None
    draw_line = None
    
    for enum, line in enumerate(lines_read):
    
        if str(line) == "% switch off outline pairs or bases by removing these lines\n":
            at_line = enum
            print(at_line)
        if str(line) == "drawoutline\n":
            draw_line = enum
            
    file_.close()
    draw_list = ["/invert false def\n", "drawreliability\n", "0.1 0.1 colorbar\n"]
    lines_read[at_line:at_line] = lines
    lines_read[draw_line + len(lines):draw_line + len(lines)] = draw_list
    
    file_ = open(file_name, "w")
    for line in lines_read:
        file_.write(line)
        
    file_.close()
    

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


    cmdline_parser = argparse.ArgumentParser('retrieve RNA seq structure')

    cmdline_parser.add_argument('-f', '--function_file',
                                default="",
                                help='input file',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-t', '--target_folder',
                                default="",
                                help='input file',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-s', '--structure_folder',
                                default="",
                                help='input file',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-a', '--summary',
                                default="",
                                help='output file summary',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-z', '--type',
                                default="",
                                help='output file summary',
                                required = True,
                                type=str)
                                
    args, unknowns = cmdline_parser.parse_known_args()
    
    

    in_path = args.target_folder
    fasta_files_to_run = glob.glob(args.structure_folder + "*")
    
    #dms_path = "/media/sven/Elements/5UTR_finished/5UTR_not_dedup_10_cov/95_10_rna_files/95_10_rna/"
    #untreated_val_path = "/media/sven/Elements/5UTR_finished/5UTR_not_dedup_10_cov/95_10_rna_files/t1/95_10_untreated_0.1/" 
    
    
    print(fasta_files_to_run)
    
    for constraint_file in fasta_files_to_run:
    
        dms_file_treated = constraint_file.split("/")[-1].split("treated")[0] + "treated.fa"
        dms_file_untreated = constraint_file.split("/")[-1].split("treated")[0] + "untreated.fa"
    
        file_ = open(constraint_file)
        lines=file_.readlines()
        id_only = lines[0]
        sequence_only = lines[1]
        constraint_only= lines[2]
        constraint_file_mod_path = "./" + constraint_file.split("/")[-1].split(".fa")[0]+ "_mod.fa"
        
        constraint_file_mod = open(constraint_file_mod_path, "w")
        constraint_file_mod.write(id_only)
        constraint_file_mod.write(sequence_only)
        constraint_file_mod.close()
        sub_string = constraint_file.split("/")[-1].split("_")[0]

        untreated_file = constraint_file.split("/")[-1].split("treated")[0] + "untreated_" + untreated_val_path.split("_")[-1].split("/")[0] +  ".fa"

        untreated_path = untreated_val_path + untreated_file


        RNAfold("RNAfold",constraint_file_mod_path, False)
        

        
        current_files = glob.glob("./*.ps")
        ind = [i for i, s in enumerate(current_files) if sub_string in s]
        sequence_name = in_path + current_files[ind[0]].split("ss")[0] +"_only_seq_"+ "ss.ps"
        os.rename(current_files[ind[0]], in_path + current_files[ind[0]].split("ss")[0] +"_only_seq_"+ "ss.ps")
                
        RNAfold("RNAfold",constraint_file, True)
        

        
        current_files = glob.glob("./*.ps")
        ind = [i for i, s in enumerate(current_files) if sub_string in s]
        treated_name = in_path + current_files[ind[0]].split("ss")[0] +"_treated_"+ "ss.ps"
        os.rename(current_files[ind[0]], in_path + current_files[ind[0]].split("ss")[0] +"_treated_"+ "ss.ps")
        color_ps_file(treated_name, dms_path + dms_file_treated.split("_constraints")[0] + "_treated.fa")
        
        RNAfold("RNAfold",untreated_path, True)
        

        
        current_files = glob.glob("./*.ps")
        ind = [i for i, s in enumerate(current_files) if sub_string in s]
        untreated_name = in_path + current_files[ind[0]].split("ss")[0] +"_untreated_"+ "ss.ps"
        os.rename(current_files[ind[0]], in_path + current_files[ind[0]].split("ss")[0] +"_untreated_"+ "ss.ps")
        font = ImageFont.truetype(args.type, size=30)
        color_ps_file(untreated_name, dms_path + dms_file_untreated.split("_constraints")[0] + "_untreated.fa")
        
        
        
        
        os. remove(constraint_file_mod_path)
        psimage1=Image.open(sequence_name)
        psimage1 = expand_white_space(psimage1)
        draw1 = ImageDraw.Draw(psimage1)
        
        l1 = draw1.textlength(constraint_file.split("/")[-1].split("_constraints")[0], font = font)
        draw1.text((psimage1.width/2 - l1/2, 5), constraint_file.split("/")[-1].split("_constraints")[0] ,(0,0,0), font = font)
        
        l1 = draw1.textlength("sequence only", font = font)
        draw1.text((psimage1.width/2 - l1/2, 40), "sequence only" ,(0,0,0), font = font)
        
        
        psimage2=Image.open(treated_name)
        psimage2 = expand_white_space(psimage2)
        draw2 = ImageDraw.Draw(psimage2)
        l2 = draw1.textlength("1lbeta-treated", font = font)
        draw2.text((psimage2.width/2 - l2/2, 30), "1lbeta-treated" ,(0,0,0), font = font)
        
        
        psimage3=Image.open(untreated_name)
        psimage3 = expand_white_space(psimage3)

        draw3 = ImageDraw.Draw(psimage3)
        l3 = draw3.textlength("untreated", font = font)
        draw3.text((psimage3.width/2 -l3/2, 30), "untreated" ,(0,0,0), font = font)
        
        dst = get_concat_v(psimage1, psimage2, psimage3)
        dst.save(in_path + constraint_file.split("/")[-1].split("treated")[0]  + ".jpg")

        os.remove(sequence_name)
        os.remove(treated_name)
        os.remove(untreated_name)
        
        log_files = glob.glob("./*.log")
        for log in log_files:
            try:
                os.remove(log)
            except:
                print("already gone")
        
        
    summary_file = open(args.summary, "w")
    summary_file.write("finished")
    summary_file.close()

    


