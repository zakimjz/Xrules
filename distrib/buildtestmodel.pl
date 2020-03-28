#!/usr/bin/perl

#usage: buildtestmodel -- -t trainDBprefix -T testDBprefix 
#                         -s supX -S supY [-n skipnewtest] 
#                         [ -d dratio -Dx xor_ratio -Dy yor_ratio]
#                         [-r runmissed -m supX -M supY] [-g skipminegetsup]
#example: -t DB/us1924/us1924 -T DB/us2430/us2430.asc -s 0.0025 -S 0.0035 -d 0.5
#for timing

$DIR="./distrib";
$XMINE = "./distrib/Xminer/xminer";

$skipnewtest=0;
$runmissed=0;
$skipmine=0;
$skipdefaultclass=0;
$msupX=0;
$msupY=0;
$dratio=0;
$or_ratio=0;
$xor_ratio=$yor_ratio = 0;
$numrules=50000;
$use_lhood_ratio=0;
$omit_vr_item=999999999; #omit vr_item is true by default
$weightmodel=2;
$skipor = 0;
$use_norm_ratio=0;

$scoringtype = 0; #avgconf
$kval = 1;

sub printusage{
    die "usage: buildtestmodel\
                -t trainDBprefix \
                -T testDBprefix \
                -s supX -S supY \
                -l use_lhood_ratio? \
                -w <weightmodel (0=prop, 1=equal, 2=inverse)> \
                [-d dratio -Dx xor_ratio -Dy yor_ratio ] \
                [-r runmissed -m supX -M supY] \
                [-n skipnewtest? ] \
                [-N use_norm_ratio? ] \
                [-C scoringtype (0=avg, 1=best, 2=bestk <k-val> )] \
                [-v omit_vitrual_item ] \
                [-x skipdefaultclass] \
                [-g skipmine? ]\n";
}

if ($#ARGV < 0){
    printusage;
}
else{
#while getopts ":d:D:gl:m:M:nrR:s:S:t:T:v:w:x" Option
    # Initial declaration.
    # a, b, c, d, e, f, and g are the options (flags) expected.
    # The : after option 'e' shows it will have an argument passed with it.
    for ($i = 0; $i <= $#ARGV; ++$i){
      switch: for ($ARGV[$i]){
          /\-C/ && do { 
              $scoringtype=$ARGV[++$i]; #def=0
              if ($scoringtype == 2){ $kval = $ARGV[++$i]; }
              last;
          }; 
          /\-d/ && do { $dratio=$ARGV[++$i]; last;}; #def=0.5
          /\-Dx/ && do { $xor_ratio=$ARGV[++$i]; last;}; #def=0.5
          /\-Dy/ && do { $yor_ratio=$ARGV[++$i]; last;}; #def=0.5
          /\-g/ && do { $skipmine=1; last;}; #def=0
          /\-l/ && do { $use_lhood_ratio=1; last;}; #def = 0
          /\-m/ && do { $msupX=$ARGV[++$i]; last;}; #reqd with -r
          /\-M/ && do { $msupY=$ARGV[++$i]; last;}; #reqd with -r
          /\-n/ && do { $skipnewtest=1; last;}; #def=0
          /\-N/ && do { $use_norm_ratio=1; last; }; #def=0
          /\-r/ && do { $runmissed=1; last;}; #def=0
          /\-R/ && do { $numrules=$ARGV[++$i]; last;}; #def=20000
          /\-s/ && do { $supX=$ARGV[++$i]; last;}; #reqd
          /\-S/ && do { $supY=$ARGV[++$i]; last;}; #reqd
          /\-t/ && do { $TRAIN=$ARGV[++$i]; last;}; #reqd
          /\-T/ && do { $TEST=$ARGV[++$i]; last;}; #reqd
          /\-v/ && do { $omit_vr_item=$ARGV[++$i]; last;}; #def is true, item=0
          /\-w/ && do { $weightmodel=$ARGV[++$i]; last;}; #def = 2; inverse
          /\-x/ && do { $skipdefaultclass=1; last;}; #def=0
          /\-help/ && do { printusage; last; };
          die "unknown argv $ARGV[$i]\n";
      }
    }
}


$TESTF="$TEST.asc";
$TRAINF="$TRAIN.asc";

if ($use_norm_ratio || !$use_lhood_ratio){
    if ($dratio == 0){$dratio=0.5;}
    if ($xor_ratio==0){$xor_ratio=0.5};
    if ($yor_ratio==0){$yor_ratio=0.5};
}
else{
    if ($dratio == 0){$dratio=1;}
    if ($xor_ratio==0){$xor_ratio=1};
    if ($yor_ratio==0){$yor_ratio=1};
}

print "RUN -t $TRAIN -T $TESTF -s $supX -S $supY -d $dratio -Dx $xor_ratio -Dy $yor_ratio -n $skipnewtest -g $skipmine -r $runmissed -m $msupX -M $msupY -v $omit_vr_item -R $numrules -l $use_lhood_ratio -w $weightmodel -nr $use_norm_ratio -C $scoringtype $kval\n";

#set out vars
$FOUT="$TRAIN.s$supX\_$supY.freq"; #output from mining X
$MODELOUT="$TRAIN.s$supX\_$supY.or$xor_ratio\_$yor_ratio.rules";
@test = split (/\//,$TEST);
$FINALOUT="$TRAIN.s$supX\_$supY.or$xor_ratio\_$yor_ratio.$test[$#test].$use_lhood_ratio.final";
$DEFOUT="$FOUT.default";
$MISSEDOUT="$FOUT.missed";

if (!$skipmine){
        $execstr = "-i $TRAINF -n 2 -s $supX $supY -o -x $omit_vr_item -r $numrules > $FOUT";
        $MINE = $XMINE; 

    #produces FOUT
    print "$MINE $execstr\n";
    `$MINE $execstr`;
}

if (!$skipnewtest){
    if (!$skipdefaultclass){
        # get default class
        $execstr="-r $FOUT -t $TRAINF -cx $xor_ratio -cy $yor_ratio -w $weightmodel -d $DEFOUT";
        if ($scoringtype){
            $execstr = join(' ', $execstr, "-s $scoringtype");
            if ($scoringtype == 2){
                $execstr = join(' ', $execstr, "$kval");
            }
        }
        if (!$use_lhood_ratio){ $execstr = join(' ', $execstr, "-C");}
        if ($use_norm_ratio){ $execstr = join(' ', $execstr, "-n");}
        $execstr = join (' ', $execstr, "> $MISSEDOUT");
        print "$DIR/newtest $execstr\n";
        `$DIR/newtest $execstr`;
    }
    #
    $execstr = "-r $FOUT -t $TESTF -cx $xor_ratio -cy $yor_ratio -w $weightmodel -D $DEFOUT";
    if ($scoringtype){
        $execstr = join(' ', $execstr, "-s $scoringtype");
        if ($scoringtype == 2){
            $execstr = join(' ', $execstr, "$kval");
        }
    }
    if (!$use_lhood_ratio){ $execstr = join(' ', $execstr, "-C");}
    if ($use_norm_ratio){ $execstr = join(' ', $execstr, "-n");}
    $execstr = join (' ', $execstr, "> $FINALOUT");
    print "$DIR/newtest $execstr\n";
    `$DIR/newtest $execstr`;
}

if ($runmissed){
    $MISSEDF="$TRAIN.missed.s$supX\_$supY.d$dratio";
    $execstr="$DEFOUT $MISSEDF";
    print "$DIR/misseddb.pl $execstr\n";
    `$DIR/misseddb.pl $execstr`;

    $execstr = "$MISSEDF";
    print "$DIR/runbin $execstr\n";
    `$DIR/runbin $execstr`;
    if ($msupX == 0) {$msupX=$supX;}
    if ($msupY == 0)  {$msupY=$supY;}

    $execstr = "-t $MISSEDF -T $MISSEDF -s $msupX -S $msupY -n -d $dratio -v $omit_vr_item -l $use_lhood_ratio";
    print "$DIR/buildtestmodel $execstr\n";
    `$DIR/buildtestmodel $execstr`;
    $MISSEDRULES="$MISSEDF.s$msupX\_$msupY.or$xor_ratio\_$yor_ratio.rules";

    $execstr = "-r $MODELOUT -t $TESTF -c $dratio -R $use_lhood_ratio -M $MISSEDRULES > $FINALOUT.missed";
    print "$DIR/newtest $execstr\n";
    `$DIR/newtest $execstr`;
}

