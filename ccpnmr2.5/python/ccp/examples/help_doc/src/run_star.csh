rm -rf ../data/nmrStar/ ../data/cnsOut/

foreach num (9 10 11 12)
  echo "\nRunning ex${num}.py"
  python ex${num}.py
  sleep 2
end
