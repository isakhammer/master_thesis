alias p3='python3.8'
#shopt
shopt -s autocd # change to named directory
shopt -s cdspell # autocorrects cd misspellings
shopt -s cmdhist # save multi-line commands in history as single line
shopt -s dotglob
shopt -s histappend # do not overwrite history
shopt -s expand_aliases # expand aliases

export EDITOR="nvim"

# # Short term alias for commands
# function ppa(){
# 	mkdir -p ~/.cache/ppa
# 	touch ~/.cache/ppa/pp
# 	echo "cd $PWD && python3 $1 || cd -" > ~/.cache/ppa/pp
# }
# alias pp="cat ~/.cache/ppa/pp && bash ~/.cache/ppa/pp "

# function ppa1(){
# 	mkdir -p ~/.cache/ppa
# 	touch ~/.cache/ppa/pp1
# 	echo "cd $PWD && python3 $1 || cd -" > ~/.cache/ppa/pp1
# }
# alias pp1="cat ~/.cache/ppa/pp1 && bash ~/.cache/ppa/pp1"

# vi mode
set -o vi
export KEYTIMEOUT=1
alias c="clear"
