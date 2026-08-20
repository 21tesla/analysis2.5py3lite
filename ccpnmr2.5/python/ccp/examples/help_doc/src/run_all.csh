rm -rf ../data/nmrStar/ ../data/xEasy/ ../data/pdb/ nmrView2XEasy/ cns2Pdb/ ../data/cnsOut/

foreach num (1 2 3 4 5 6 7 8 9 10 11 12)
  echo "\nRunning ex${num}.py"
  python ex${num}.py
  sleep 2
end
