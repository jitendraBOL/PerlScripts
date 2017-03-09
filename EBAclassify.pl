#!/usr/local/bin/perl
		
########################################
#
#  Perl Program to classify the species and convert them into Evolutionary Breakpoint Analyser format- Classify 
#  It required NCBI taxonomy database ( taxdump) to classify the species in EBA format. 
#  The location of taxdump folder ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

#  Author: Jitendra Narayan < jnarayan81@gmail.com>
#  Copyright (c) 2014 by Jitendra. All rights reserved.
#  You may distribute this script under the same terms as perl itself
#

use strict;

$|++; ## flush perl’s print buffer

my (%names, %nodes, @array2d, %finalHash);

printUsage ();

## Input species names.
print "Enter your species name [Comma seperated] -- First name should be your reference species name : \n";
my $speciesString=<STDIN>;
chomp ($speciesString); 
my @species = split(/,/, $speciesString);
s{^\s+|\s+$}{}g foreach (@species); # Removing leading and trailing whitespace from array strings.
foreach (@species) { $_ =~ s/\s+/_/g; }   ## add unerscore to space.

#######

my $taxonomy_directory='taxdump';
  my ($files,$dirs)=getDirectoryFiles($taxonomy_directory);
  process_file($_) for @$files;    # This subrutine create a hash of names and nodes file.

# my %reverseNames = map {$names{lc($_)} => lc($_)} keys(%names);   # Reverse the names hash and search the inpit species name
my %reverseNames;
while (my($key, $value) = each %names) { $reverseNames{lc($value)}=$key; }

my @inName=split(/:/, $reverseNames {lc($species[0])}); 
my $inputSpsId = $inName[0]; 
print "\nYour reference species Id and Name is : $inputSpsId\t$species[0]\n";

# shift(@species);  ## No need to Remove the first elements which is Reference species becuase it is being used for phylogenetic classification.

my $num=0; my %species_id;
print "\n\tFollowings are the Ids and name of your studied species\t\n";
foreach my $sp (@species) {
	my @array; 
	my @k=findkey ($sp, \%names );  
        s{^\s+|\s+$}{}g foreach @k; # Removing leading and trailing whitespace from array strings.
	 
	# if (keys %{{ map {$_, 1} @k }} != 1) { @k[0]=senseMultipleHits(\@k,$sp);} 
	LINE: if (scalar(@k) > 1) { @k[0]=senseMultipleHits(\@k,$sp);}       
	elsif (!@k) { ($sp,@k)=senseNoHits(\@k,$sp); goto LINE;}
	else { print "\t@k\t$sp\n";}  # It search the "species names" in name file and extract its key.

 
	my @val = split(/:/, $k[0]);
	my $val = $val[0];  push @array, $val;  # Leaf of the tree ids are stored here, later we will add nodes in this array.
	%species_id; $species_id{$val}=$sp;   # stored the hash for future use ... only studied species id 

		while ($val > 1) {
		my $parent = $nodes{$val};  # print "$parent\t";
        	#my $name= $names {$parent}; print "$name\n";  ### need to correct it print the different name bcz of duplication
		$val=$parent;
		push @array, $parent;
                }

	my @new_array=reverse @array;

	for (my $aa=0; $aa<= $#new_array; $aa++) { $array2d[$aa][$num]=$new_array[$aa]; }    ### Store the data in 2d array

undef @array; $num++;
}

#print2d(@array2d);  # To print the 2d array of all species classification.
my %finalHash=classify(\@array2d, \%reverseNames);  ## the main subruitine to classify the species.
#printhash(\%names);   # To print the hash

my @finalHashValues = values %finalHash;
my %finalHash2 = merge(\@finalHashValues, \%finalHash); 
my $counter=0;
foreach my $spsKey (keys %finalHash2) {       ## print the final result form finalHash array
    if ($counter == 0 ) { print "\n==========================================\nlineage=\n";} 
    print "$spsKey=$finalHash2{$spsKey}\n";
$counter++;
}

#### All subrutines here ####

#----------------------------------------------
sub senseMultipleHits {
	my ($k_ref, $sp)=@_;
	my  @k=@$k_ref;
 	print "\t-----FATAL Warning-----\n\tWe have more than one hits in database and all of them are not the same\nPlease choose any one of them by Number[#]\n\t<#>\t<Ids>\t<Name>\n";
	my %multipleId; 
	for(my $num=0; $num<=$#k; $num++) { print "\t$num\t$k[$num]\t$sp\n"; $multipleId{$num}=$k[$num];} 
	my $opted=<STDIN>; chomp($opted); 
	$k[0]=$multipleId{$opted}; 
	print "You have opted $opted=$k[0]\t$sp\n";
return @k[0];
 ## Check the array for similarities of all hits and if they are different from one another then flash a message. 
}

#----------------------------------------------
sub senseNoHits {
  my ($k_ref, $sp)=@_;
  my  @k=@$k_ref;

  # sensorium to make sense of any errors
  print "\t----- FATAL Error -----\n\tNo hits were found for $sp in our dataset !!!\n";
  print "Lets try with a new name [ Scientific Name suggested ]: ";
  my $newName = <STDIN>; chomp($newName); $newName=trim($newName); $newName =~ s/\s+/_/g;  ## replace the space with underscore ...
  my @spsId=findkey ($newName, \%names );

return ($newName,@spsId);
}

#----------------------------------------------------------------------
# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

#----------------------------------------------------------------------
sub getDirectoryFiles{          # It get the directory files and return it
     my $taxdir = shift;

     opendir(my $dh, $taxdir) || die "can't opendir $taxdir : $!";
     my @entries = grep {!( /^\.$/ || /^\.\.$/)} readdir($dh);
     @entries =  map { "$taxdir/$_" } @entries; #change to absolute paths
     closedir $dh;

     my @files =  grep( -f $_ , @entries);
     my @dirs = grep(-d $_, @entries);
     return (\@files,\@dirs);     ## return as a reference 


}

#----------------------------------------------------------------------

sub process_file{ # This is your custom subroutine to perform on each file    
    my $f = shift;     
    my ($val, $nam) = check_file($f);   
	if ($val == 1 and ($nam eq "names"))
		{  # print "processing file $f\n";
		%names=file2hash ($f, $nam);
		#return @array;
		}
         elsif ($val == 1 and $nam eq "nodes")
		{ # print "processing file $f\n";
		%nodes=file2hash ($f, $nam);
		}
        
	}

#-----------------------------------------------------------------------

sub check_file {
    use File::Basename;
    my $filepath = shift; # print $file;
    my $file = basename($filepath);
    my @ff= split /\./, $file;
    if ($ff[0] eq "names" || "nodes" ) 
	{ return 1, $ff[0]; }
}

#-------------------------------------------------------------------------
sub file2hash {
	my ($infile, $n) = @_;
	my %hash;
	open FILE, $infile or die $!;
	while (<FILE>)
		{
   		chomp;  # s/^\s*(.*)\s*$/$1/;
		my @tmp_array= split /\t/ , $_;   
		s{^\s+|\s+$}{}g foreach @tmp_array; # Removing leading and trailing whitespace from array strings.
   		my ($key, $val) = split /\t\|\t/;
   		# if (($infile eq "names.dmp") and  ($tmp_array[6] ne "scientific name")) {next;} # print "$tmp_array[6]\n";      ## if we want to enter only specific lines.
   		if($n eq "names") { $key="$key:$.";}    # I make it unique by adding the line number and split later...
		$val =~ s/\s+/_/g;  ## replace the space with underscore ...
		$hash{$key} = $val;  
		# print "$n\t$key\t$val\n";
        } 
	close FILE;
return %hash;      
}

#-------------------------------------------------------------------------

sub printhash {
  	my %hash=%{$_[0]};
	foreach my $key (sort keys %hash) {
     	print "$key : $hash{$key}\n";
	}
}

#-------------------------------------------------------------------------
sub findkey {
        my ($species, $hash) =@_;
	my  %hash=%$hash;  my @all_keys;
	foreach my $key (keys %hash) {
     	if ($hash{$key} =~ m/^$species$/i) { push @all_keys, $key};
	}
s{^\s+|\s+$}{}g foreach @all_keys; # Removing leading and trailing whitespace from array strings.
return @all_keys;
undef @all_keys;
} 

#-------------------------------------------------------------------------
sub findId {
        my ($id, $hash) =@_;
	my  %hash=%$hash;  my @all_values;
	foreach my $key (keys %hash) {
	my @newValue= split(/:/, $hash{$key});       ## What if we have more than two hits for a key !!!!!
     	if ($newValue[0]==$id) { push @all_values, $key};
	}
s{^\s+|\s+$}{}g foreach @all_values; # Removing leading and trailing whitespace from array strings.
my $all_values=join(",",uniq(@all_values));
return $all_values;
undef @all_values;
} 

#-----------------------------------------------------------------------
sub print2d {
	my @array_2d=@_;
	for(my $i = 0; $i <= $#array_2d; $i++){
	   for(my $j = 0; $j <= $#{$array_2d[0]} ; $j++){
	      print "$array_2d[$i][$j]\t";
	   }
	   print "\n";
	}
}

#---------------------------------------------------------------------
sub classify {
      my ($array2d_ref, $reverseNames_ref)=@_; 
      my  @array2d=@$array2d_ref;
      my  %reverseNames=%$reverseNames_ref;
      my @column;   my @all_leaf; my @row; my @uarray;  my @dup; my @all_g;  my $flag=0;

      for(my $i = 0; $i <= $#array2d; $i++){
	   for(my $j = 0; $j <= $#{$array2d[0]} ; $j++){
	      push @row, $array2d[$i][$j];

		#for (@{$array2d[$i]}) { print "$_\t";}  ## it store all the values of $i rows.
		@uarray=uniq(@{$array2d[$i]}); 
		my %seen; @dup = map { 1==$seen{$_}++ ? $_ : () } @{$array2d[$i]};   ## to search only duplicated values in an array;
		@dup = grep {$_} @dup;   ## delete the blank or 0 in an array;
		# print "@dup\t---@uarray\t====\n";

                my @column = map {$$_[$j]} @array2d;   
		@column = grep {$_} @column; # delete the blank array at the end.
		push (@all_leaf, $column[-1]);
			}
                	@all_leaf=uniq(@all_leaf);		
			#if (keys %{{ map {$_, 1} @dup }} == 1)
			if(scalar(@dup) == scalar(@uarray))
				{
				my @breaks_in= class2eba (\@all_leaf, \@all_leaf, $inputSpsId); 
				# print "@dup: @breaks_in\n";
				my @breakInNames =  returnName (\%species_id, \@breaks_in); # print @breakInNames;
				if ($flag != 1) {  my $brkin=join(",",uniq (@breakInNames)); $finalHash{$inputSpsId}= $brkin;}     
				$flag=1;
				}
 			else
				{
				my %all_clustered = clustered (\@dup,\@array2d);
				foreach my $group (keys %all_clustered) {
    					#print "The members of $group are\n";
					#no strict 'refs'; print "$group\n";
					my @uall=uniq(@{$all_clustered{$group}});
   					my @breaks_in= class2eba (\@uall, \@all_leaf, $inputSpsId);
					my @breakIn =  returnName (\%species_id, \@breaks_in);
					my $breaks_in=join(",",uniq (@breakIn)); ## print all sub breakpoints species
					#my $groupName = findId ($group, \%reverseNames);   ### group is not identified in species_id because we have ids of only studied species.  Can use to extract name !!!!
					$finalHash{$group}= $breaks_in;
   					#print "$group:@breaks_in\n";
					push @all_g, $group;
    					}
				}
	undef @all_leaf; #print "$_\n";
	}

my $all_class = join(", ",@all_g);
$finalHash{'classification'}="lineage,$all_class,$inputSpsId";
return %finalHash;
undef @all_g;
}
 #closedir(DIR);
 
##---------------------------------------------------------------------
sub merge {       # It merge the key from hash
        my ($values_array,$hash) =@_;
	my @values_array=@$values_array;
	my  %hash=%$hash;  
       	my %result; my @all_keys;
      
        for my $values (@values_array) {
		for my $key (keys %hash) {

        	my $accu = $hash{$key};
           	if ($values eq $accu) { push ( @all_keys, $key); }
	
             }
	my $keystring=join(",",@all_keys);

        undef @all_keys;
	$result{$keystring} = $values;

    }
return %result;
}

##------------------------------------------------------------------
sub returnName
{
my ($speciesId_ref,$breakIn_ref)= @_;
my %speciesId = %$speciesId_ref;
my @breakInName;
my @breakIn = @$breakIn_ref;
	foreach my $ids (@breakIn) { 
		my $idsName =  $speciesId { $ids };
		push (@breakInName, $idsName);
	#	print "$idsName\t";
		}
return @breakInName;
undef @breakInName;
}

#--------------------------------------------------------------------
	
sub uniq { my %seen; return grep { !$seen{$_}++ } @_; }

#---------------------------------------------------------------------
# Checks if a provided element exists in the provided list
# Usage: isInList <needle element> <haystack list>
# Returns: 0/1
sub isInList {
    my $needle = shift;
    my @haystack = @_;
    foreach my $hay (@haystack) {
        if ( $needle eq $hay ) {
            return 1;
        }
    }
    return 0;
}

	 
#-----------------------------------------------------------------------
sub clustered {
   my ($dup_ref,$array2d_ref)= @_;
   my @dup = @$dup_ref;
   my @array2d = @$array2d_ref;

   my @group; my %all_group;

	foreach my $d(@dup) { 

		for(my $i = 0; $i <= $#array2d; $i++){
	   		for(my $j = 0; $j <= $#{$array2d[0]} ; $j++){
	      			
               		 	my @column = map {$$_[$j]} @array2d;   
				@column = grep {$_} @column; # delete the blank array at the end.
				my $val = isInList ($d, @column);
				if ($val == 1) {
				push (@group, $column[-1]);  
						}
					}
				}
	$all_group{$d}=[@group]; 
	undef @group;			 
		}
return %all_group;
}

#----------------------------------------------------------------------------

sub class2eba  {
	my ($leaf_ref,$all_leaf_ref,$reference)= @_;
   	my @leaf = @$leaf_ref;
	my @all_leaf = @$all_leaf_ref;
	chomp($reference);
	my $val = isInList($reference, @leaf);
	 #print "$val\n";
		if($val != 1) { return @leaf;}
		else { 
	
			my %in_leaf = map {$_ => 1} @leaf;
			my @diff  = grep {not $in_leaf{$_}} @all_leaf;	
			       if ($#diff != -1) { return @diff; }
				else {
					my @new_leaf = grep { $_ != $reference } @all_leaf;
					return @new_leaf;
				}
		    }

	}
	
#}


