#!/bin/bash

# Initialiser le compteur d'échecs à 0
count_errors=0

# Nombre de fois que tu veux exécuter le programme
iterations=100000

# Boucle pour exécuter le programme plusieurs fois
for (( i=1; i<=iterations; i++ ))
do
    # Exécuter le programme (remplace "ton_programme" par le nom de ton programme)
    make test_index_multi > "mtest_$i.out" 2> "mtest_$i.out"

    # Vérifier si le programme a échoué (code de retour différent de 0)
    if [ $? -ne 0 ]; then
        echo "Échec détecté lors de l'exécution $i"
        count_errors=$((count_errors + 1))  # Incrémenter le compteur d'échecs
        cat "mtest_$i.out" | tail -n 10000 > "mtest_$i.final"
        cat "mtest_$i.final" | tail -n 100
        rm "mtest_$i.out"
        # exit 0
    else
        echo -ne "Pas probleme $i\r"
        rm "mtest_$i.out"
    fi
done

# Afficher le nombre total d'échecs
echo "Nombre total d'échecs : $count_errors"
