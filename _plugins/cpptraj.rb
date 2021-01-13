# -*- coding: utf-8 -*- #
# frozen_string_literal: true

## Use the Bash one as basis
# module Rouge
#   module Lexers

## WILL NOT WORK on GitHub Pages
# This "hook" is executed right before the site's pages are rendered
# https://stackoverflow.com/questions/61814907/how-to-specify-a-custom-language-parser-alias-for-rouge-in-jekyll-3

Jekyll::Hooks.register :site, :pre_render do |site|
  puts "Adding cpptraj aliases..."
  require "rouge"

  module Rouge
    module Lexers

    class Cpptraj < RegexLexer
      title "Cpptraj"
      desc 'Cpptraj, the AMBER program for processing trajectories and data'
      tag 'cpptraj'
      aliases 'traj'
      filenames '*.in'

      KEYWORDS = %w(
        trajin trajout autoimage center reference rms parm
      ).join('|')

      BUILTINS = %w(
        lastframe ref out onlyframes origin
      ).join('|')

      state :basic do
        rule %r/#.*$/, Comment

        rule %r/\b(#{KEYWORDS})\s*\b/, Keyword
        rule %r/\bcase\b/, Keyword, :case

        rule %r/\b(#{BUILTINS})\s*\b(?!(\.|-))/, Name::Builtin
        rule %r/[.](?=\s)/, Name::Builtin

        rule %r/(\b\w+)(=)/ do
          groups Name::Variable, Operator
        end

        rule %r/[\[\]{}()!=>]/, Operator
        rule %r/&&|\|\|/, Operator

        # here-string
        rule %r/<<</, Operator

        rule %r/(<<-?)(\s*)(['"]?)(\\?)(\w+)(\3)/ do |m|
          groups Operator, Text, Str::Heredoc, Str::Heredoc, Name::Constant, Str::Heredoc
          @heredocstr = Regexp.escape(m[5])
          push :heredoc
        end
      end

      state :heredoc do
        rule %r/\n/, Str::Heredoc, :heredoc_nl
        rule %r/[^$\n\\]+/, Str::Heredoc
        mixin :interp
        rule %r/[$]/, Str::Heredoc
      end

      state :heredoc_nl do
        rule %r/\s*(\w+)\s*\n/ do |m|
          if m[1] == @heredocstr
            token Name::Constant
            pop! 2
          else
            token Str::Heredoc
          end
        end

        rule(//) { pop! }
      end


      state :double_quotes do
        # NB: "abc$" is literally the string abc$.
        # Here we prevent :interp from interpreting $" as a variable.
        rule %r/(?:\$#?)?"/, Str::Double, :pop!
        mixin :interp
        rule %r/[^"`\\$]+/, Str::Double
      end

      state :ansi_string do
        rule %r/\\./, Str::Escape
        rule %r/[^\\']+/, Str::Single
        mixin :single_quotes
      end

      state :single_quotes do
        rule %r/'/, Str::Single, :pop!
        rule %r/[^']+/, Str::Single
      end

      state :data do
        rule %r/\s+/, Text
        rule %r/\\./, Str::Escape
        rule %r/\$?"/, Str::Double, :double_quotes
        rule %r/\$'/, Str::Single, :ansi_string

        # single quotes are much easier than double quotes - we can
        # literally just scan until the next single quote.
        # POSIX: Enclosing characters in single-quotes ( '' )
        # shall preserve the literal value of each character within the
        # single-quotes. A single-quote cannot occur within single-quotes.
        rule %r/'/, Str::Single, :single_quotes

        rule %r/\*/, Keyword

        rule %r/;/, Punctuation

        rule %r/--?[\w-]+/, Name::Tag
        rule %r/[^=\*\s{}()$"'`;\\<]+/, Text
        rule %r/\d+(?= |\Z)/, Num
        rule %r/</, Text
        mixin :interp
      end

      state :curly do
        rule %r/}/, Keyword, :pop!
        rule %r/:-/, Keyword
        rule %r/[a-zA-Z0-9_]+/, Name::Variable
        rule %r/[^}:"`'$]+/, Punctuation
        mixin :root
      end

      # the state inside $(...)
      state :paren_interp do
        rule %r/\)/, Str::Interpol, :pop!
        rule %r/\(/, Operator, :paren_inner
        mixin :root
      end

      # used to balance parentheses inside interpolation
      state :paren_inner do
        rule %r/\(/, Operator, :push
        rule %r/\)/, Operator, :pop!
        mixin :root
      end

      state :math do
        rule %r/\)\)/, Keyword, :pop!
        rule %r([-+*/%^|&!]|\*\*|\|\|), Operator
        rule %r/\d+(#\w+)?/, Num
        mixin :root
      end

      state :case do
        rule %r/\besac\b/, Keyword, :pop!
        rule %r/\|/, Punctuation
        rule %r/\)/, Punctuation, :case_stanza
        mixin :root
      end

      state :case_stanza do
        rule %r/;;/, Punctuation, :pop!
        mixin :root
      end

      state :backticks do
        rule %r/`/, Str::Backtick, :pop!
        mixin :root
      end

      state :interp do
        rule %r/\\$/, Str::Escape # line continuation
        rule %r/\\./, Str::Escape
        rule %r/\$\(\(/, Keyword, :math
        rule %r/\$\(/, Str::Interpol, :paren_interp
        rule %r/\${#?/, Keyword, :curly
        rule %r/`/, Str::Backtick, :backticks
        rule %r/\$#?(\w+|.)/, Name::Variable
        rule %r/\$[*@]/, Name::Variable
      end

      state :root do
        mixin :basic
        mixin :data
      end
    end
    end
end
end
