
alias p3='python3.8'
#shopt
shopt -s autocd # change to named directory
shopt -s cdspell # autocorrects cd misspellings
shopt -s cmdhist # save multi-line commands in history as single line
shopt -s dotglob
shopt -s histappend # do not overwrite history
shopt -s expand_aliases # expand aliases

export EDITOR="nvim"
alias juplab="jupyter-lab --ip 0.0.0.0 --port 8888 --allow-root"

# vi mode
set -o vi
export KEYTIMEOUT=1
alias c="clear"
