#!/bin/bash

cat > .input.aq <<EOF
molecule
	coords cartesian
	atom O  0.00000000  0.00000000  0.11726921
	atom H  0.75698224  0.00000000 -0.46907685
	atom H -0.75698224  0.00000000 -0.46907685
	basis
		basis_set cc-pVDZ

scf
	convergence 1e-9
	max_iterations 100
	damping 0
	diis
		order 6
		start 8

cc
	diis
		order 5
	convergence 1e-9
	max_iterations 100

cholesky
	delta 1e-10
	cond_max 1e3
EOF

echo
echo 'Testing AO RHF CCSD + Lambda-CCSD'
echo

(../bin/bench-ao-ccsd-lambda .input.aq | tee .out) || echo 'Executable ../bin/bench-ao-ccsd-lambda not found' && exit 1

cat .out | \
    perl -0777 -ne 'print $1 if /SCF(.*)Orbital/s' | \
    perl -0777 -ne 'print "$1\n" while m/([0-9]+ +-?[0-9]+\.[0-9]+ +[0-9]+\.[0-9]+[eE][\+-][0-9]+)/gms' | \
    tail -n 1 |
{
    read it en conv
    conv=`echo $conv | sed 's/[eE]/\*10\^/g'`
    
    echo -n "SCF: "
    
    if [ `echo "$conv > 10^-9" | bc -l` -eq 1 ]; then
        echo "failed - not converged ($conv > 10^-9)"
        exit 1
    fi
    
    diff=`echo "($en)-(-74.55012645669238)" | bc -l`
    if [ `echo "$diff >  10^-9" | bc -l` -eq 1 -o `echo "$diff< -10^-9" | bc -l` -eq 1 ]; then
        echo "failed - energy differs ($diff)"
        exit 1
    fi
    
    echo "passed"
} || exit 1

cat .out | perl -0777 -ne 'print $1 if /MP2 Energy: +(-?[0-9]+\.[0-9]+)/' |
{
    read en
    
    echo -n "MP2: "
    
    diff=`echo "($en)-(-0.171348679568)" | bc -l`
    if [ `echo "$diff >  10^-9" | bc -l` -eq 1 -o `echo "$diff< -10^-9" | bc -l` -eq 1 ]; then
        echo "failed - energy differs ($diff)"
        exit 1
    fi
    
    echo "passed"
} || exit 1

cat .out | \
    perl -0777 -ne 'print $1 if /CCSD(.*)Lambda/s' | \
    perl -0777 -ne 'print "$1\n" while m/([0-9]+ +-?[0-9]+\.[0-9]+ +[0-9]+\.[0-9]+[eE][\+-][0-9]+)/gms' | \
    tail -n 1 |
{
    read it en conv
    conv=`echo $conv | sed 's/[eE]/\*10\^/g'`
    
    echo -n "CCSD: "
    
    if [ `echo "$conv > 10^-9" | bc -l` -eq 1 ]; then
        echo "failed - not converged ($conv > 10^-9)"
        exit 1
    fi
    
    diff=`echo "($en)-(-0.180145524753)" | bc -l`
    if [ `echo "$diff >  10^-9" | bc -l` -eq 1 -o `echo "$diff< -10^-9" | bc -l` -eq 1 ]; then
        echo "failed - energy differs ($diff)"
        exit 1
    fi
    
    echo "passed"
} || exit 1

cat .out | \
    perl -0777 -ne 'print $1 if /Lambda(.*)Final/s' | \
    perl -0777 -ne 'print "$1\n" while m/([0-9]+ +-?[0-9]+\.[0-9]+ +[0-9]+\.[0-9]+[eE][\+-][0-9]+)/gms' | \
    tail -n 1 |
{
    read it en conv
    conv=`echo $conv | sed 's/[eE]/\*10\^/g'`
    
    echo -n "Lambda: "
    
    if [ `echo "$conv > 10^-9" | bc -l` -eq 1 ]; then
        echo "failed - not converged ($conv > 10^-9)"
        exit 1
    fi
    
    diff=`echo "($en)-(-0.180145524753)" | bc -l`
    if [ `echo "$diff >  10^-9" | bc -l` -eq 1 -o `echo "$diff< -10^-9" | bc -l` -eq 1 ]; then
        echo "failed - energy differs ($diff)"
        exit 1
    fi
    
    echo "passed"
} || exit 1

rm -f .input.aq .out
exit 0
