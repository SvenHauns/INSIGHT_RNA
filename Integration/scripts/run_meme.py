import subprocess as sp
import glob
import os
import argparse






def run_meme(input_sequences, input_background_file, output_html):


	cmd = "meme {targets} -dna -oc {meme_out} -neg {background} -nmotifs 15 -minw 6 -maxw 15 -revcomp -objfun de"
	cmd = cmd.format(targets=input_sequences, background=input_background_file, meme_out = output_html)
	
	log_file = "meme.log"
	print(cmd)
	with open(log_file, 'w') as lf:
		sp.call(cmd.split(), stdout=lf)
        
    
	return
    
    

    
    
if __name__ == "__main__":


    cmdline_parser = argparse.ArgumentParser('retrieve RNA seq structure')

    cmdline_parser.add_argument('-f', '--input_sequences',
                                default="",
                                help='input file',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-t', '--input_background_file',
                                default="",
                                help='background file',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-y', '--output_html',
                                default="",
                                help='output html',
                                required = True,
                                type=str)

                                
    args, unknowns = cmdline_parser.parse_known_args()
    run_meme(args.input_sequences, args.input_background_file, args.output_html)

