_micronota_completion() {
    COMPREPLY=( $( env COMP_WORDS="${COMP_WORDS[*]}" \
                   COMP_CWORD=$COMP_CWORD \
                   _MICRONOTA_COMPLETE=complete $1 ) )
    return 0
}

complete -F _micronota_completion -o default micronota;
