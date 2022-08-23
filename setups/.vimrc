" show line numbers
set number

"Highlight matching search patterns
set hlsearch
" Enable incremental search
set incsearch
" Include matching uppercase words with lowercase search term
set ignorecase
" Include only uppercase words with uppercase search term
set smartcase

" Display options
set showmode
set showcmd


" Enable 256 colors
set t_Co=256

" Color scheme
" colo desert

" Enable mouse usage (all modes)
set mouse=a

filetype plugin indent on
" show existing tab with 4 spaces width
set tabstop=4

" when indenting with '>', use 4 spaces width
set shiftwidth=4

" On pressing tab, insert 4 spaces
set expandtab

" Turn off vim bell
set visualbell

"au FocusGained,BufEnter * :checktime
