#!/bin/bash

# Variables
first_command="make test_index_mono"  # Remplace par ta commande initiale
second_command="make test_index_multi"  # Remplace par ta deuxième commande EN MULTI
output_file="resultat_initial.txt"
temp_file="resultat_temp.txt"
iterations=100  # Nombre de répétitions pour la deuxième commande

# Exécuter la première commande et sauvegarder le résultat dans un fichier
echo "Exécution de la première commande..."
$first_command | tail -n +3 > "$output_file"
echo "Résultat initial sauvegardé dans $output_file."

# Boucle pour exécuter plusieurs fois la deuxième commande
for ((i = 1; i <= iterations; i++)); do
    # echo "Exécution de la commande secondaire, tentative $i..."
    $second_command | tail -n +3 > "$temp_file" 2> "$temp_file"

    # Comparer les fichiers
    if cmp -s <(grep -v "^##" "$output_file") <(grep -v "^##" "$temp_file"); then
        echo "Tentative $i : Les fichiers sont identiques."
    else
        echo "Tentative $i : Les fichiers sont différents.";
        diff <(grep -v "^##" "$output_file") <(grep -v "^##" "$temp_file");
        # exit 0;
    fi
done

# Nettoyage
rm "$temp_file"
