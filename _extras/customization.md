---
layout: page
title: Customizing Your Workshop's Website
permalink: /customization/index.html
---
## Configuration File `_config.yml`

You should edit the `_config.yml` configuration file in the root directory of your workshop to
configure some site-wide variables and make the site function correctly:

* `carpentry` - to tell us which carpentry workshop this is, possible values are ("swc" for Software
  Carpentry workshops, "dc" for Data Carpentry workshops, "lc" for Library Carpentry Workshops, or
  "cp" for general Carpentries events such as instructor trainings).
* `curriculum` - for Data Carpentry, which one of the curriculum is being taught. Possible values
  are: `dc-ecology`, `dc-genomics`, `dc-socsci`, `dc-geospatial`.
* `flavor` - `r` or `python` depending on which lessons are being taught at the workshop (currently
  only for Data Carpentry workshops)
* `title` - overall title for the workshop. If set (i.e., different from "Workshop Title" or empty),
  it will appear in the "jumbotron" (the gray box at the top of the page). This variable is also
  used for the title of the extra pages. More information about extra pages are [available in the
  README](https://github.com/carpentries/workshop-template#creating-extra-pages).

For example, if the URL for the repository is `https://github.com/gvwilson/2015-07-01-miskatonic`,
the URL for the website will be `http://gvwilson.github.io/2015-07-01-miskatonic`.

You should not need to modify any of the other values in `_config.yml`.

## Home Page (`index.md`): data in the YAML header

Your workshop's home page lives in `index.md`,
which must define the values below in its header.
If your workshop is taught online, see the
[following section](#for-online-workshops) for customization
options.

*   `layout` must be `workshop`.

*   `venue` is the short name of the institution or group hosting the
    workshop, like "Euphoric State University".  It should *not*
    include the address or other details, since this value is
    displayed in a table on websites (e.g.,
    <https://carpentries.org/upcoming_workshops/>). See section
    below for value to use for online workshops.

*   `address` is the workshop's address (including details like the
    room number). The address should be all on one line.
    See section below for value to use for online workshops.

*   `country` must be a two-letter ISO-3166 code for the country in
    which the workshop is going to take place, such as "fr" (for
    France) or "nz" (for New Zealand) - see [Wikipedia](https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2#Officially_assigned_code_elements)
    for a complete list. See section below for value to use for 
    online workshops.

*   `language` is the language that will be used in the workshop.
    It must be an [ISO 639-1 code](https://en.wikipedia.org/wiki/List_of_ISO_639-1_codes).
    Note that two-letter codes mean different things for countries
    and languages: "ar" is Arabic when used for a language, but
    Argentina when used for a country.

*   `latitude` and `longitude` are the latitude and longitude of the workshop
    site (so we can put a pin on our map). You can use
    [this site](https://getlatlong.net/) to find these values. 
    See section below for value to use for online workshops.

*  `humandate` is the human-friendly start and end date for the
    workshop.  Please use three-letter month names and abbreviations
    (e.g., `Jul` instead of `July`), since these values are displayed
    in a table on our websites.  (Strictly speaking this information
    is redundant, since we require a machine-readable `startdate` and
    `enddate`, but reliably translating those into human-readable
    dates is an interesting challenge...)

*   `humantime` is the human-friendly start and end time for each day of
    the workshop, e.g., "09:00 am - 4:00 pm" or "09:00-16:00".  (We
    recognize that we ought to allow different start or end times on
    different days, but going down that path leads eventually to
    embedding iCal date/time specifications in our headers, which in
    turn leads to madness...)

*   `startdate` is the workshop's starting date in YYYY-MM-DD format,
    such as `2015-07-01`.  You must use four digits for the year and
    two each for the month and day.

*   `enddate` is the workshop's ending date in the same format.  If your
    workshop is only one day long, the `enddate` field should be deleted.
    If your workshop has a more complicated schedule (e.g., a half day a
    week for four weeks), please delete the `enddate` field and only tell
    us its start date.

*   `instructor` is a comma-separated list of instructor names.  The
    list must be enclosed in square brackets, and each name must be in
    double quotes, as in `["Alan Turing","Grace Hopper"]`.  Do not
    include other information (such as the word "instructor") in these
    values.

*   `helper` is a comma-separated list of helper names formatted in the
    same way as the instructor names.  If there are no helpers, use an
    empty list `[]`.

*   `contact` is the contact email address to use for your workshop.
    If you do not provide a contact email address, your website will
    display the address for the workshop coordinators (who probably
    won't be able to answer questions about the specific details of
    your workshop).

The header may optionally define the following:

*   `collaborative_notes` is the URL for the Etherpad for your workshop.
    If you are not using an Etherpad, you can delete this line. You can
    create a carpentries etherpad [here](https://pad.carpentries.org/).

*   `eventbrite` is the multi-digit Eventbrite registration key.  If you
    are using Eventbrite, the Carpentries Regional Coordinators will
    give this to you.  If you are using something else, you may delete
    this line.  Note: this value must be given as a string in double
    quotes, rather than as a number.

### For online workshops

If the workshop is online, follow the same instructions as above with the
following modifications:

* `venue`: Use the name of the institution that organizes the workshop and do
  not include a mention that it is an online workshop.
* `address`: If you can safely share the URL for the videoconferencing, you may
  list it here (it must start with `http` or `https`); if you cannot or prefer
  to not share the videoconferencing information, use the value `online`.
* `country`: Please use the country associated with the host institution for the
  workshop.
* `latitude` and `longitude`: if it makes sense, use the coordinates for the
  host institution. If it does not, use `0` for both the latitude and the
  longitude.


## Home Page: Schedule and Syllabus

You should edit the sections titled `Schedule` and `Syllabus`
so that they show what you're actually planning to teach and when.  These
files are located in the appropriate workshop folder (`dc`, `lc` or `swc`)
inside the `_includes` folder.

## Home Page: Setup

You may have to edit the `setup.html` located
in the `dc`, `lc` or `swc` folders
inside the `_includes` folder.
Make sure you delete the sections for the
software and data
you will **not** be using in your workshop,
so that learners don't spend time installing
software they don't need.


## Updating the repository

If for some reason,
such as the installation instructions having become disconnected
with the current lesson material,
you need to get changes from this repository
into the clone of it with your workshop page,
please follow the steps bellow:

1.  Add the workshop-template repository as upstream:

        $ git remote add upstream https://github.com/carpentries/workshop-template.git

2.  Fetch the data from upstream repository (also know as the workshop-template
    repository):

        $ git pull upstream

4.  Address possible merge conflicts, and

        $ git commit -a

5.  Push the changes to your repository on GitHub:

        $ git push origin gh-pages