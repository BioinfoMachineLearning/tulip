echo "Script started for TULIP"
a=(1rdi 1cx9 1akv 1str 2d06 )
for i in "${a[@]}"; do
  echo "Running with protein: $i"
  python3 complex_generator.py --protein_name="$i" --config=../configs/configs.yml &
done
wait
echo "All proteins are done"


