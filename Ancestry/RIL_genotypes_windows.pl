#!/usr/bin/perl -w
#use strict;

my $OutputFile = 'FR_RIL_genotypes_winZI1k.txt';

my $WinFilePrefix = 'windows_ZI1k_Chr';

my @chrs = ('X','2L','2R','3L','3R');

#my @RILs = ("E101R","E102R","E103R","E104R","E105R","E106R","E107R","E108R","E109R","E110R","E111R","E112R","E113R","E115R","E116R","E117R","E118R","E119R","E121R","E122R","E123R","E124R","E126R","E127R","E128R","E12R","E130R","E131R","E135R","E136R","E137R","E139R","E13R","E141R","E142R","E143R","E144R","E145R","E147R","E148R","E149R","E14R","E150R","E152R","E153R","E154R","E155R","E156R","E157R","E158R","E159R","E160R","E163R","E164R","E165R","E168R","E169R","E170R","E173R","E174R","E180R","E181R","E182R","E183R","E184R","E185R","E186R","E187R","E188R","E189R","E18R","E190R","E191R","E192R","E193R","E194R","E195R","E197R","E198R","E199R","E19R","E1R","E200R","E201R","E202R","E203R","E204R","E205R","E206R","E207R","E208R","E209R","E20R","E212R","E213R","E214R","E218R","E219R","E221R","E222R","E223R","E224R","E225R","E226R","E227R","E228R","E229R","E22R","E230R","E232R","E234R","E236R","E237R","E238R","E23R","E241R","E242R","E245R","E246R","E248R","E249R","E24R","E250R","E251R","E253R","E254R","E255R","E256R","E258R","E259R","E260R","E261R","E262R","E264R","E265R","E266R","E267R","E269R","E26R","E271R","E272R","E273R","E274R","E275R","E277R","E278R","E279R","E280R","E282R","E283R","E284R","E285R","E286R","E287R","E288R","E289R","E290R","E291R","E292R","E293R","E294R","E295R","E297R","E29R","E2R","E300R","E303R","E305R","E306R","E307R","E308R","E309R","E30R","E310R","E311R","E313R","E315R","E316R","E317R","E319R","E320R","E321R","E322R","E326R","E327R","E328R","E329R","E32R","E330R","E331R","E332R","E333R","E334R","E336R","E337R","E338R","E339R","E33R","E340R","E341R","E342R","E343R","E345R","E346R","E347R","E348R","E34R","E351R","E352R","E353R","E354R","E356R","E358R","E359R","E360R","E361R","E362R","E364R","E365R","E366R","E367R","E368R","E369R","E36R","E371R","E372R","E374R","E375R","E377R","E378R","E379R","E37R","E380R","E381R","E382R","E384R","E385R","E386R","E38R","E39R","E3R","E40R","E41R","E42R","E43R","E44R","E45R","E46R","E47R","E48R","E4R","E50R","E51R","E52R","E53R","E54R","E55R","E57R","E58R","E60R","E61R","E62R","E63R","E65R","E66R","E69R","E6R","E70R","E71R","E72R","E73R","E74R","E76R","E78R","E7R","E80R","E81R","E82R","E83R","E84R","E85R","E86R","E87R","E88R","E89R","E8R","E90R","E91R","E92R","E93R","E94R","E95R","E96R","E97R","E98R","E9R");

my @RILs = ("F100R","F101R","F102R","F103R","F104R","F105R","F106R","F107R","F108R","F109R","F110R","F111R","F113R","F114R","F115R","F117R","F118R","F119R","F11R","F120R","F121R","F122R","F123R","F124R","F125R","F126R","F127R","F128R","F129R","F12R","F130R","F131R","F132R","F133R","F134R","F135R","F136R","F137R","F13R","F140R","F141R","F142R","F143R","F144R","F145R","F146R","F147R","F148R","F149R","F14R","F151R","F152R","F153R","F154R","F155R","F157R","F158R","F159R","F15R","F160R","F161R","F162R","F163R","F165R","F166R","F167R","F168R","F169R","F16R","F170R","F172R","F173R","F175R","F176R","F177R","F178R","F179R","F17R","F181R","F182R","F183R","F184R","F186R","F187R","F188R","F18R","F190R","F191R","F194R","F195R","F197R","F199R","F19R","F200R","F201R","F203R","F205R","F206R","F207R","F208R","F209R","F20R","F210R","F212R","F213R","F214R","F215R","F216R","F217R","F218R","F219R","F21R","F220R","F221R","F222R","F226R","F227R","F228R","F230R","F231R","F232R","F234R","F235R","F236R","F237R","F238R","F23R","F240R","F241R","F242R","F243R","F244R","F246R","F247R","F248R","F249R","F24R","F251R","F252R","F253R","F257R","F258R","F259R","F25R","F260R","F262R","F263R","F264R","F265R","F266R","F267R","F270R","F271R","F273R","F275R","F276R","F278R","F279R","F27R","F281R","F282R","F283R","F284R","F285R","F286R","F287R","F28R","F291R","F292R","F293R","F294R","F297R","F298R","F299R","F300R","F301R","F302R","F303R","F304R","F305R","F306R","F307R","F308R","F309R","F310R","F311R","F312R","F313R","F314R","F315R","F316R","F317R","F318R","F319R","F31R","F320R","F321R","F323R","F324R","F325R","F326R","F327R","F329R","F32R","F330R","F331R","F332R","F333R","F334R","F336R","F33R","F341R","F342R","F343R","F345R","F346R","F349R","F34R","F350R","F351R","F352R","F353R","F354R","F358R","F359R","F360R","F361R","F363R","F366R","F368R","F369R","F36R","F371R","F372R","F374R","F375R","F378R","F381R","F382R","F383R","F385R","F386R","F38R","F391R","F393R","F395R","F398R","F3R","F401R","F402R","F403R","F404R","F405R","F407R","F408R","F409R","F40R","F410R","F41R","F42R","F43R","F44R","F45R","F46R","F47R","F48R","F49R","F50R","F51R","F53R","F54R","F55R","F56R","F57R","F58R","F59R","F5R","F60R","F62R","F63R","F64R","F65R","F69R","F71R","F73R","F74R","F76R","F78R","F79R","F81R","F83R","F84R","F85R","F86R","F88R","F8R","F90R","F91R","F92R","F93R","F94R","F95R","F96R","F97R","F98R","F99R");

my @handles = @RILs;
for ($i = 0; $i < @handles; $i++){
  $handles[$i] =~ s/1/A/;
  $handles[$i] =~ s/2/B/;
  $handles[$i] =~ s/3/C/;
  $handles[$i] =~ s/4/D/;
  $handles[$i] =~ s/5/E/;
  $handles[$i] =~ s/6/F/;
  $handles[$i] =~ s/7/G/;
  $handles[$i] =~ s/8/H/;
  $handles[$i] =~ s/9/I/;
  $handles[$i] =~ s/0/J/;
}

my $c = 0;
my $f = 0;
my $i = 0;
my $j = 0;
my $r = 0;
my $s = 0;
my $w = 0;
my $chr = '';
my $file = '';
my $site = 0;
my $Win2Sums = 0;
my $Win1Sums = 0;
my $Win0Sums = 0;

my @line = ();
my @WinStarts = ();
my @WinStops = ();
my @InputAoA = ();
my @ChrGenoAoA = ();

#open output file
open O, ">$OutputFile";
print O "chr\tWinStart\tWinStop";
for ($r = 0; $r < @RILs; $r++){
  print O "\t$RILs[$r]";
}
print O "\n";

#for each chr arm
for ($c = 0; $c < @chrs; $c++){
  $chr = $chrs[$c];
  @WinStarts = ();
  @WinStops = ();
  @ChrGenoAoA = ();
  @line = ();
  for ($w = 0; $w < @WinStops; $w++){
    push @ChrGenoAoA, [ @line ];
  }
  
#open window file
  $file = $WinFilePrefix . $chr . '.txt';
  open W, "<$file" or die;
  while (<W>){
    chomp;
    last if m/^$/;
    @line = split;
    push @WinStarts, $line[0];
    push @WinStops, $line[1];
  }
  close W;
  $w = 0;
  
#open and analyze each RIL file
  for ($r = 0; $r < @RILs; $r++){
    @InputAoA = ();
    $file = $RILs[$r] . '.posterior';
    open $handles[$r], "<./$chr/$file" or die "can not open ./$chr/$file\n";
    $file = $handles[$r];
    scalar (<$file>);
    while (<$file>){
      chomp;
      last if m/^$/;
      @line = split;
      push @InputAoA, [ @line ];
    }
    close $file;
    $site = 0;
    for ($w = 0; $w < @WinStops; $w++){
      $Win2Sums = 0;
      $Win1Sums = 0;
      $Win0Sums = 0;
      for ($s = $site; $s < @InputAoA; $s++){
	$site = $s;
	last if ($InputAoA[$s][1] > $WinStops[$w]);
	$Win2Sums += $InputAoA[$s][2];
	$Win1Sums += $InputAoA[$s][3];
	$Win0Sums += $InputAoA[$s][4];
      }
      if (($Win2Sums + $Win1Sums + $Win0Sums) == 0){
	push @{$ChrGenoAoA[$w]}, -999;
	next;
      }
      if ($Win2Sums > $Win1Sums){
	if ($Win2Sums > $Win0Sums){
	  push @{$ChrGenoAoA[$w]}, 2;
	}
	else{
	  push @{$ChrGenoAoA[$w]}, 0;
	}
      }
      else{
	if ($Win1Sums > $Win0Sums){
	  push @{$ChrGenoAoA[$w]}, 1;
	}
	else{
	  push @{$ChrGenoAoA[$w]}, 0;
	}
      }
###
#      print "$RILs[$r] $site $ChrGenoAoA[$w][-1] $InputAoA[$s][2] $InputAoA[$s][3] $InputAoA[$s][4]\n";
#      last;
###      
    }
    print "Finished chr $chr, RIL $RILs[$r]\n";
  }

###
#  die;
  $j = @ChrGenoAoA;
  print "$j x ";
  $j =  @{$ChrGenoAoA[0]};
  print "$j\n";
#  for ($w = 0; $w < @ChrGenoAoA; $w++){
#    for ($r = 0; $r < @{$ChrGenoAoA[0]}; $r++
###
    
#add chr genotypes to output
  for ($w = 0; $w < @WinStops; $w++){
    print O "$chr\t$WinStarts[$w]\t$WinStops[$w]";
    for ($r = 0; $r < @{$ChrGenoAoA[0]}; $r++){
      print O "\t$ChrGenoAoA[$w][$r]"
    }
    print O "\n";
  }
}
close O;


    
