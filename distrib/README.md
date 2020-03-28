
1) run makeall


2) Type buildtestmodel.pl to see the options:

usage: buildtestmodel

                -t trainDBprefix
                -T testDBprefix
                -s supX -S supY
                -l use_lhood_ratio?
                -w <weightmodel (0=prop, 1=equal, 2=inverse)>
                [-d dratio -Dx xor_ratio -Dy yor_ratio ]
                [-r runmissed -m supX -M supY]
                [-n skipnewtest? ]
                [-N use_norm_ratio? ]
                [-C scoringtype (0=avg, 1=best, 2=bestk <k-val> )]
                [-v omit_vitrual_item ]
                [-x skipdefaultclass]
                [-g skipmine? ]

Example run:

 buildtestmodel.pl -t synth.train -T synth.test -s 0.2 -S 0.2
 
	where -t gives the training file (prefix only; .asc is assumed)
	      -T gives the test file (prefix only; .asc is assumed)
	      -s support for class X (works for 2 class problems only)
	      -S support for class Y (works for 2 class problems only)
		 NOTE: buildtestmodel works for 2 classes, 
		  but Xminer can mine any number of classes
        now you must choose the strength measure:
	   default is confidence, if you do not give any other option
	   specify -l for likelihood ratio
		   -N for weighted confidence

	now -v will allow Xminer to omit the root node in the XML tree
	    -v X will exculde X from the output

	YOU MUST USE -v 0 for the cslogs dataset (the root item is a 
	             fake item to join a forest into a tree)
	       

	You can modify the ratio to prune rules (the min rule strength)
	   for (weighted) conf the default is 0.5 
	   while for likelihood is it 1.0 
           you can give any value by specifying
	   -d -Dx -Dy
	   e.g.: -d 2 -Dx 2 -Dy 2

3) run for example 

	buildtestmodel.pl -t synth.train -T synth.test -s 0.2 -S 0.2 -l

   the final output is in the file:

    synth.train.s0.25_0.25.or1_1.synth.train.1.final
	
   There are 4 sets of values output:

   1st entry is for the partial classifier, without the default 
       class determination  

   2nd entry is for proportional cost model

   3rd entry is for equal cost model
   
   4th entry is for inverse cost model

   e.g.:

        DEFAULT PREDICTIONS, LABEL= classX 2
        class    True    False
        X       9       0
        Y       0       0
        TRUE/FALSE INFO FOR NPREDICT+DEFAULT
        class   True    False   Accuracy        Coverage
        X       55      8       0.873016        0.873016
        Y       58      8       0.878788        0.878788
        TOTAL   63      66      0.875835        1

This means that for this cost model the default class was X
It then gives the number of rules that do not have any rule, 
so must be classified as the default class (9 in this case)

It then gives the confusion matrix and the accuracy and 
coverage for each class and for the whole classifier.
