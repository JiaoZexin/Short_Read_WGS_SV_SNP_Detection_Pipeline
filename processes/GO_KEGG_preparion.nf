process GO_KEGG_preparion {
    tag "go_kegg"

    input:
    path pep_fasta

    output:
    path "12.GO_KEGG/*"

    script:
    """
    mkdir -p 12.GO_KEGG
    cd 12.GO_KEGG

    # Step 1: Download go basic 
    wget -q http://purl.obolibrary.org/obo/go/go-basic.obo

    # Step 2: Using HuaSmallTools 
    git clone https://github.com/Hua-CM/HuaSmallTools.git
    python HuaSmallTools/parse/parse_go_obofile.py -i go-basic.obo -o go.tb

    # Step 3: emapper
    emapper.py -m diamond -i ../${pep_fasta} -o output_pep --output_dir . --cpu 4

    # Step 4: Generate GO/KEGG background
    python HuaSmallTools/parse/parse_eggNOG.py \
      -i output_pep.emapper.annotations \
      -g go.tb \
      -o output_pep.GOKEGG.DB

    echo "GOannotation.tsv and KOannotation.tsv ready as enrichment background." > complete.log
    """
}