process Fimpute_imputation {
    tag "fimpute"

    input:
    path control_file

    output:
    path "13.fimpute/*"

    script:
    """
    mkdir -p 13.fimpute
    cd 13.fimpute

    start_time=\$(date +%s)
    echo "[`date`] Running FImpute3 imputation" > fimpute.log

    # Modify based on your enviroment
    FImpute3 control_file_final.txt >> fimpute.log 2>&1

    end_time=\$(date +%s)
    echo "FImpute3 finished in \$((end_time - start_time)) seconds" >> fimpute.log
    """
}