using CSV, DataFrames, ProgressMeter, BioSequences, BioAlignments, Random, Statistics

function complement(nt) 
    if nt == "A" 
        return "T"
    elseif nt == "T" 
        return "A"
    elseif nt == "C" 
        return "G"
    elseif nt == "G" 
        return "C"
    else
        return nt  # Pour gérer les cas N ou autres
    end
end

function reverse_complement(seq)
    return join(reverse([complement(x) for x in split(seq,"")]))
end

function read_fasta_ref(filename)
    mutants = []
    sequences = []
    io = open(filename, "r")
    while !eof(io)
        header = readline(io)
        sequence = readline(io)
        # Nettoyer le header en enlevant le ">" initial s'il existe
        header = startswith(header, ">") ? header[2:end] : header
        push!(mutants, header)
        push!(sequences, sequence)
    end
    close(io)
    return DataFrame(mutant=mutants, sequence=sequences)
end

function read_fastq(filename)
    sequences = []
    quality_scores = []
    headers = []
    
    io = open(filename, "r")
    while !eof(io)
        header = readline(io)
        sequence = readline(io)
        separator = readline(io)
        quality = readline(io)
        
        # Nettoyer le header en enlevant le "@" initial et tout ce qui suit le premier espace
        header = startswith(header, "@") ? header[2:end] : header
        header = split(header, " ")[1]
        
        # Calculer le score de qualité
        quality_score = sum([(Int(x)-33) for x in quality])
        
        push!(headers, header)
        push!(sequences, sequence)
        push!(quality_scores, quality_score)
    end
    
    close(io)
    return sequences, quality_scores, headers
end

function compare_seq(seq1, seq2, decalage=0)
    # Vérifier que les séquences ne sont pas vides après le décalage
    if length(seq1) <= decalage
        return [], 1.0
    end
    
    seq1_trim = seq1[decalage+1:end]
    seq2_trim = seq2[1:min(length(seq1_trim), length(seq2))]
    
    # Ajuster la longueur des séquences pour qu'elles soient égales
    min_length = min(length(seq1_trim), length(seq2_trim))
    seq1_trim = seq1_trim[1:min_length]
    seq2_trim = seq2_trim[1:min_length]
    
    missmatch = []
    for i in 1:min_length
        if seq1_trim[i] != seq2_trim[i]
            push!(missmatch, i)
        end
    end
    
    taux_erreur = length(missmatch)/min_length
    return missmatch, taux_erreur
end

function find_closest_mutant(seq::AbstractString, df_ref::DataFrame, decalage::Int64=0)
    best_score = 0.0
    best_mutant = ""
    best_missmatch = []
    
    for i in 1:nrow(df_ref)
        ref_seq = df_ref.sequence[i]
        missmatch, taux_erreur = compare_seq(seq, ref_seq, decalage)
        score = 1.0 - taux_erreur
        
        if score > best_score
            best_score = score
            best_mutant = df_ref.mutant[i]
            best_missmatch = missmatch
        end
    end
    
    return best_mutant, best_score, best_missmatch
end

function compare_read_1(seq_exp_1::Vector{Any}, df_ref::DataFrame, n::Int64, decalage::Int64=0, quality_scores::Vector{Any}=[], headers::Vector{Any}=[])
    if n > length(seq_exp_1)
        n = length(seq_exp_1)
        println("n trop grand, on prend n = ", n)
    end

    # Créer le DataFrame pour les résultats
    df_results = DataFrame(
        Header = String[],
        Mutant_1 = String[],
        Score_1 = Float64[],
        Sequence_1 = String[],
        Quality_1 = Int[]
    )

    # Sélectionner n indices aléatoires
    id = randperm(length(seq_exp_1))[1:n]
    
    # Barre de progression
    p = Progress(n, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50)
    
    for i in id
        seq_1_t = seq_exp_1[i]
        quality_score = !isempty(quality_scores) ? quality_scores[i] : 0
        header = !isempty(headers) ? headers[i] : string(i)
        
        mutant1, score1, missmatch = find_closest_mutant(seq_1_t, df_ref, decalage)
        
        push!(df_results, (
            header,
            mutant1,
            score1,
            seq_1_t,
            quality_score
        ))
        
        next!(p)
    end

    return df_results
end

function main()
    # Paramètres
    n = 100000
    base_dir = dirname(dirname(@__FILE__))  # Chemin absolu vers APE007_analysis
    output_dir = joinpath(base_dir, "data", "phe_analysis")
    
    # Vérification et création des dossiers
    if !isdir(output_dir)
        @info "Création du dossier de sortie : $output_dir"
        mkpath(output_dir)
    end
    
    # Fichiers à analyser avec chemins absolus
    data_dir = joinpath(base_dir, "data")
    init_file = joinpath(data_dir, "APE007-31-08-23-init_S9_L001_R1_001.fastq")
    phe_file = joinpath(data_dir, "APE007-31-08-23-phe_S10_L001_R1_001.fastq")
    ref_file = joinpath(data_dir, "S1Aref.fasta")
    
    # Vérification de l'existence des fichiers
    for (desc, file) in [("initial", init_file), ("phe", phe_file), ("référence", ref_file)]
        if !isfile(file)
            error("Le fichier $desc n'existe pas : $file")
        end
    end
    
    # Lecture des fichiers
    println("Lecture des fichiers de séquences...")
    println("Lecture du fichier de référence...")
    seq_ref = read_fasta_ref(ref_file)
    println("Nombre de séquences de référence : ", nrow(seq_ref))
    
    println("\nLecture des séquences initiales...")
    seq_init, quality_init, headers_init = read_fastq(init_file)
    println("Nombre de séquences initiales : ", length(seq_init))
    
    println("\nLecture des séquences phe...")
    seq_phe, quality_phe, headers_phe = read_fastq(phe_file)
    println("Nombre de séquences phe : ", length(seq_phe))
    
    # Analyse des séquences
    println("\nAnalyse des séquences initiales...")
    df_init = compare_read_1(seq_init, seq_ref, n, 24, quality_init, headers_init)
    init_output = joinpath(output_dir, "resultat_init.csv")
    CSV.write(init_output, df_init)
    println("Résultats initiaux sauvegardés dans : $init_output")
    
    println("\nAnalyse des séquences phe...")
    df_phe = compare_read_1(seq_phe, seq_ref, n, 25, quality_phe, headers_phe)
    phe_output = joinpath(output_dir, "resultat_phe.csv")
    CSV.write(phe_output, df_phe)
    println("Résultats phe sauvegardés dans : $phe_output")
    
    # Afficher un résumé des résultats
    println("\nRésumé des résultats :")
    println("Nombre de mutants uniques (init) : ", length(unique(df_init.Mutant_1)))
    println("Nombre de mutants uniques (phe) : ", length(unique(df_phe.Mutant_1)))
    println("Score moyen (init) : ", mean(df_init.Score_1))
    println("Score moyen (phe) : ", mean(df_phe.Score_1))
    
    println("\nAnalyse terminée. Les résultats sont disponibles dans le dossier : $output_dir")
end

# Exécution du programme
main() 