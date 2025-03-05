using CSV, DataFrames, Statistics

function extraire_phenotype(mutant_name)
    # Le format est : "organisme,classe,ordre,famille,cold/hot,phenotype"
    try
        champs = split(mutant_name, ",")
        if length(champs) >= 6
            return champs[6]  # Le phénotype est le dernier champ
        end
        return "inconnu"
    catch
        return "inconnu"
    end
end

function calculer_enrichissement(df_init::DataFrame, df_phe::DataFrame)
    # Créer un dictionnaire pour compter les occurrences
    counts_init = Dict()
    counts_phe = Dict()
    phenotypes = Dict()
    
    # Compter les occurrences dans init
    for mutant in df_init.Mutant_1
        counts_init[mutant] = get(counts_init, mutant, 0) + 1
        phenotypes[mutant] = extraire_phenotype(mutant)
    end
    
    # Compter les occurrences dans phe
    for mutant in df_phe.Mutant_1
        counts_phe[mutant] = get(counts_phe, mutant, 0) + 1
        phenotypes[mutant] = extraire_phenotype(mutant)
    end
    
    # Créer le DataFrame avec le format demandé
    df_enrichissement = DataFrame(
        mutant = String[],
        phenotype = String[],
        count_init = Int[],
        count_sort = Int[],
        enrichment = Float64[]
    )
    
    # Calculer l'enrichissement pour chaque mutant
    mutants_uniques = unique(vcat(collect(keys(counts_init)), collect(keys(counts_phe))))
    
    for mutant in mutants_uniques
        count_i = get(counts_init, mutant, 0)
        count_s = get(counts_phe, mutant, 0)
        
        # Calculer l'enrichissement (log2 fold change)
        # Ajouter un pseudo-count de 1 pour éviter les divisions par zéro
        enrichissement = log2((count_s + 1) / (count_i + 1))
        
        push!(df_enrichissement, (
            mutant,
            phenotypes[mutant],
            count_i,
            count_s,
            enrichissement
        ))
    end
    
    # Trier par enrichissement décroissant
    sort!(df_enrichissement, :enrichment, rev=true)
    
    return df_enrichissement
end

function main()
    # Définition des chemins
    base_dir = dirname(dirname(@__FILE__))  # Chemin absolu vers APE007_analysis
    data_dir = joinpath(base_dir, "data", "phe_analysis")
    
    # Vérification du dossier de données
    if !isdir(data_dir)
        error("Le dossier d'analyse n'existe pas : $data_dir")
    end
    
    # Chemins des fichiers
    init_file = joinpath(data_dir, "resultat_init.csv")
    phe_file = joinpath(data_dir, "resultat_phe.csv")
    
    # Vérification des fichiers d'entrée
    for (desc, file) in [("initial", init_file), ("phe", phe_file)]
        if !isfile(file)
            error("Le fichier de résultats $desc n'existe pas : $file")
        end
    end
    
    # Charger les données
    println("Chargement des données...")
    df_init = CSV.read(init_file, DataFrame)
    df_phe = CSV.read(phe_file, DataFrame)
    
    println("Nombre de séquences dans init : ", nrow(df_init))
    println("Nombre de séquences dans phe : ", nrow(df_phe))
    
    # Calculer l'enrichissement
    println("\nCalcul de l'enrichissement...")
    df_enrichissement = calculer_enrichissement(df_init, df_phe)
    
    # Afficher un résumé des résultats
    println("\nRésumé des résultats :")
    println("Nombre total de mutants : ", nrow(df_enrichissement))
    println("Nombre de phénotypes uniques : ", length(unique(df_enrichissement.phenotype)))
    println("Phénotypes trouvés : ", join(unique(df_enrichissement.phenotype), ", "))
    
    # Sauvegarder les résultats
    println("\nSauvegarde des résultats...")
    output_file = joinpath(data_dir, "enrichissement.csv")
    CSV.write(output_file, df_enrichissement)
    println("Résultats sauvegardés dans : $output_file")
    
    println("\nAnalyse terminée. Les résultats sont disponibles dans le dossier : $data_dir")
end

# Exécution du programme
main() 