"""Build the saved standalone note; optionally update the user's existing editor file."""
from pathlib import Path
import argparse
import shutil

ROOT=Path(__file__).resolve().parents[1]
P=argparse.ArgumentParser()
P.add_argument('--editor',type=Path)
args=P.parse_args()
bit=(ROOT/'research/independent/bit-improvement/construction.tex').read_text()
bit=bit.split(r'\section{Assembly choices with a linear phase guard}')[0]
guard=(ROOT/'research/independent/guard-improvement/semantic_guard.tex').read_text()
guard=guard[guard.index(r'\section{The distinction'):guard.index(r'\paragraph{Reproduction and attribution.}')]
preamble=r'''\documentclass[11pt]{article}
\usepackage[margin=1in]{geometry}
\usepackage[T1]{fontenc}
\usepackage{amsmath,amssymb,amsthm,booktabs,hyperref}
\hypersetup{colorlinks=true,urlcolor=blue,linkcolor=blue}
\newtheorem{lemma}{Lemma}
\newtheorem{proposition}{Proposition}
\newtheorem{theorem}{Theorem}
\newcommand{\ii}{\mathrm i}
\newcommand{\Z}{\mathbb Z}
\title{Conditional integer-multiplication bounds\\with a linear precision guard}
\author{A research extension prepared with OpenAI Codex}
\date{October 8, 2026}
\begin{document}
\maketitle
'''
content=preamble+(ROOT/'notes/latest-intro.tex').read_text()+'\n'+guard+'\n'+bit+'\n'+(ROOT/'notes/latest-assembly.tex').read_text()+'\n\\end{document}\n'
target=ROOT/'notes/complex-reuse-note.tex'
target.write_text(content)
if args.editor:
 if not args.editor.is_file():raise SystemExit('The existing editor file must already exist.')
 args.editor.write_text(content)
print('Built standalone note:',target)
