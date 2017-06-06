#!/usr/bin/perl -w

# This program uses a list of enzymatic reactions and prints in the standard
# output all systems of Differential-Algebraic Equations (DAEs) that may be 
# obtained through different simplifications based on all possible combinations
# of replacement of Ordinary Differential Equations (ODEs) by mass conservation
# algebraic relations.
#
#    Copyright (C) 2017 Marcelo S. Reis.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

# If you use this software, please cite:
#
# Marcelo S. Reis, Vincent Noël, Matheus H.S. Dias, Layra L. Albuquerque,
# Amanda S. Guimarães, Lulu Wu, Junior Barrera, and Hugo A. Armelin.
#
#       "An interdisciplinary approach for designing kinetic models
#                  of the Ras/MAPK signaling pathway." 
#
# Kinase Signaling Networks, Methods in Molecular Biology, vol. 1636, chap. 28.
# doi: 10.1007/978-1-4939-7154-1_28.
#

use strict;


# Set this flag as 1 if it is supposed to use Michaelis-Menten equation for all
# reactions and 0 otherwise.
#
my $MICHAELIS_MENTEN_SIMPLIFICATION_FOR_ALL_REACTIONS = 1;


# Set this flag as 1 if it is supposed to add a feedback of pp-ERK on Raf* and
# 0 otherwise.
#
my $ppERK_FEEDBACK = 1; 


# Each index of this array corresponds to one reaction of Eq. 1 and 2.1--2.9,
# precisely in that order. for instance, $reactions[0] is linked to the 
# enzymatic reaction of Eq. 1:
#
#                RasGTP + Raf k-1 <-> k1 RasGTP-Raf -> k1cat RasGTP + Raf*
# 
my @reactions = ( ['RasGTP' , 'Raf'    , 'Raf*'],      
                  ['Pase1'  , 'Raf*'   , 'Raf'],
                  ['Raf*'   , 'MEK'    , 'p-MEK'],
                  ['Raf*'   , 'p-MEK'  , 'pp-MEK'],
                  ['Pase2'  , 'p-MEK'  , 'MEK'],
                  ['Pase2'  , 'pp-MEK' , 'p-MEK'],
                  ['pp-MEK' , 'ERK'    , 'p-ERK'],
                  ['pp-MEK' , 'p-ERK'  , 'pp-ERK'],
                  ['Pase3'  , 'p-ERK'  , 'ERK'],
                  ['Pase3'  , 'pp-ERK' , 'p-ERK']                  
                  );

if ($ppERK_FEEDBACK)
{
  $reactions[10] = ['pp-ERK' , 'Raf*' , 'Raf']; 
}

# Algebraic Equations (AEs) that define the mass conservation relations.
# Raf0, MEK0 and ERK0 denote the total quantity of Raf, MEK and ERK,
# respectively.
#                 
my %AE;

if ($MICHAELIS_MENTEN_SIMPLIFICATION_FOR_ALL_REACTIONS)
{
  # When we adopt the Quasi-steady-state (QSS) assumption, we can
  # measure either:
  #
  # [MEK0] and [p-MEK] + [pp-MEK]
  #
  # or 
  # 
  # [ERK0] and [p-ERK] + [pp-ERK]
  #
  %AE = (
   'Raf0' => ['Raf', 'Raf*'],
   'MEK0' => ['MEK', 'p-MEK',  'pp-MEK'],
   'ERK0' => ['ERK', 'p-ERK', 'pp-ERK'] 
        );
}
else
{
  # Without the QSS assumption, in our time-course Western blot assays,
  # we assume we can measure either: 
  #
  # [MEK0] and [p-MEK] + [pp-MEK] + [pp-MEK-ERK] + [pp-MEK-p-ERK] + [Raf*-p-MEK] 
  #
  # or 
  #
  # [ERK0] and [p-ERK] + [pp-ERK] + [pp-MEK-p-ERK] + [pp-ERK-Raf*] 
  #
  %AE = (
   'Raf0' => ['Raf', 'Raf*', 'Raf*-MEK', 'Raf*-p-MEK', 'pp-ERK-Raf*'],
   'MEK0' => ['MEK', 'p-MEK', 'pp-MEK', 'pp-MEK-ERK', 'pp-MEK-p-ERK',
              'Raf*-MEK', 'Raf*-p-MEK'],
   'ERK0' => ['ERK', 'p-ERK', 'pp-ERK', 'pp-MEK-ERK', 'pp-MEK-p-ERK'] 
        );
        
  if ($ppERK_FEEDBACK)
  {
    $AE{'ERK0'} = ['ERK', 'p-ERK', 'pp-ERK', 'pp-ERK-Raf*', 
                   'pp-MEK-ERK', 'pp-MEK-p-ERK'];
  }
}


# System of Ordinary Differential Equations (ODEs):
#   
my %ODE;

my $id = 0;


# Construct the original ODE system:
#
foreach my $reaction (@reactions)
{
  $id++;

  my $enzyme    = $reaction->[0];
  my $substrate = $reaction->[1];
  my $product   = $reaction->[2];
  
  if (($MICHAELIS_MENTEN_SIMPLIFICATION_FOR_ALL_REACTIONS)
     || 
     (($enzyme eq 'RasGTP')  || ($enzyme eq 'Pase1') || ($enzyme eq 'Pase2') 
      || ($enzyme eq 'Pase3')))
  {
    # For enzymes that are not dependent variables, we consider them
    # as constants and use the Michaelis-Menten equation.
    #
    $ODE{$substrate} .= " - kcat$id" . "[$enzyme][$substrate]/K$id" 
                        . "m+[$substrate]"; 

    $ODE{$product} .= " + kcat$id" . "[$enzyme][$substrate]/K$id"
                        . "m+[$substrate]"; 

  }
  else
  {
    $ODE{$enzyme} .= " - k$id" . "[$enzyme][$substrate] + k-$id" 
                    . "[$enzyme-$substrate] - k$id" . "cat[$enzyme-$substrate]"; 

    $ODE{$substrate} .= " - k$id" . "[$enzyme][$substrate] + k-$id" 
                        . "[$enzyme-$substrate]"; 

    $ODE{"$enzyme-$substrate"} .= " + k$id" . "[$enzyme][$substrate] - k-$id" 
                    . "[$enzyme-$substrate] - k$id" . "cat[$enzyme-$substrate]"; 

    $ODE{$product} .= " + k$id" . "cat[$enzyme-$substrate]"; 
  }
}


# Printing the original ODE system:
#
my $system_size = 0;
foreach my $species (keys %ODE)
{
  printf "d[$species]/dt =%s\n\n", $ODE{$species};
  my @terms = split " ", $ODE{$species};
  $system_size += (@terms / 2);
}

print "Size of the original ODE system: $system_size right-side terms.\n\n\n";


# Now we print all the possible simplified DAE systems:
#
my $system_id = 0;

foreach my $Raf_species (@{$AE{'Raf0'}})
{
  my $Raf_AE = "[$Raf_species] = [Raf0]";
  
  foreach my $right_side_species (@{$AE{'Raf0'}})
  {
    $right_side_species ne $Raf_species 
      and $Raf_AE .= " - [$right_side_species]";
  }

  foreach my $MEK_species (@{$AE{'MEK0'}})
  {
    my $MEK_AE = "[$MEK_species] = [MEK0]";
  
    foreach my $right_side_species (@{$AE{'MEK0'}})
    {
      $right_side_species ne $MEK_species 
        and $MEK_AE .=  " - [$right_side_species]";
    }
 
    foreach my $ERK_species (@{$AE{'ERK0'}})
    {
      my $ERK_AE = "[$ERK_species] = [ERK0]";
  
      foreach my $right_side_species (@{$AE{'ERK0'}})
      {
        $right_side_species ne $ERK_species 
          and $ERK_AE .=  " - [$right_side_species]";
      }

      my $system_size = 0;
      foreach my $species (keys %ODE)
      {
        if ($species eq $Raf_species)
        {
          print $Raf_AE . "\n\n";
        }
        elsif ($species eq $MEK_species)
        {
          print $MEK_AE . "\n\n";
        }
        elsif ($species eq $ERK_species)
        {
          print $ERK_AE . "\n\n";
        }
        else
        {
          printf "d[$species]/dt =%s\n\n", $ODE{$species};
          my @terms = split " ", $ODE{$species};
          $system_size += (@terms / 2);
        }
      }
      $system_id++;
      print "Size of DAE system $system_id: $system_size" .
            " right-side terms.\n\n\n";

    } # Raf + MEK + ERK

  } # Raf + MEK

} # Raf

# End of program.
#
exit 0;
 

