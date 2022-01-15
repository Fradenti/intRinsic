## Resubmission

This is a resubmission. In this version we have:

* Shorten the DESCRIPTION title to comply with the 67 characters limit.

* Added the references in the description field of our DESCRIPTION. The adopted
format is Author (Year, <doi:...>). This modification implied the identification of possibly misspelled words (last names of authors, "et", and "al"). Please, see below for more details.
  
* Added a \value{} section to all .Rd files, addressing the 11 missing Rd-tags.

* Changed the \dontrun{} wrapper, previously used for lengthy examples, to \donttest{}.
  
## R CMD check results


Running `devtools::check(args = c('--as-cran','--no-manual'))` locally produces 
no ERRORs, WARNINGs, nor NOTEs.  

Running `devtools::check_rhub()` produces no ERRORs, no WARNINGs, and the 
following NOTEs.  

- On most of the platforms, we obtained

* checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Francesco Denti <francescodenti.personal@gmail.com>’
  New submission

  Possibly misspelled words in DESCRIPTION:
    al (17:14, 18:14, 19:16, 20:25)
    Denti (18:5)
    et (17:11, 18:11, 19:13, 20:22)
    Facco (17:5)

Please note that the misspelled words detected are part of the references that we were asked to add.
  
- Only on `Platform: Windows Server 2022, R-devel, 64 bit`, we obtained

* checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

We have not been able to replicate this NOTE, and we do not understand why
it is produced. However, running `devtools::check_win_devel()` and
`devtools::check_win_release()` only produce

* checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Francesco Denti <francescodenti.personal@gmail.com>’
  New submission
  
  Possibly misspelled words in DESCRIPTION:
   al (17:14, 18:14, 19:16, 20:25)
   Denti (18:5)
   et (17:11, 18:11, 19:13, 20:22)
   Facco (17:5)

Finally, this package, in its current state, also passes all the standard 
checks performed via GitHub action.

## Downstread dependencies

There are currently no downstream dependencies for this package.
