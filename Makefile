all: website

# directory containing my org files as well as my assets files
SRC_DIR ?= src
# directory where I will but the files for my website (HTML + assets)
DST_DIR ?= _site

# list all files in src
# if you want to exclude .org files use the exclude from the find command
SRC_RAW_FILES := $(shell find $(SRC_DIR) -type f)
# generate all file that should be copied in the site
# For my site, I want to publish my source files along the HTML files
DST_RAW_FILES   := $(patsubst $(SRC_DIR)/%,$(DST_DIR)/%,$(SRC_RAW_FILES))
ALL             += $(DST_RAW_FILES)

# COPY EVERYTHING 
$(DST_DIR)/% : $(SRC_DIR)/%; mkdir -p "$(dir $@)"; cp "$<" "$@"

# MD -> HTML
EXT := .md
# all source file we'll pass to pandoc
SRC_PANDOC_FILES ?= $(shell find $(SRC_DIR) -type f -name "*$(EXT)")
# all destination files we expect (replace the extension by .html)
DST_PANDOC_FILES ?= $(subst $(EXT),.html, \
                        $(subst $(SRC_DIR),$(DST_DIR), \
                            $(SRC_PANDOC_FILES)))
ALL              += $(DST_PANDOC_FILES)

# use a template (you should use one)
# TEMPLATE ?= templates/post.html
# URL of the CSS put yours
CSS = /css/bb.css
# The pandoc command to run to generate an html out of a source file
PANDOC := pandoc \
            --to html5 \
			--mathjax \
            -c $(CSS) \
            --standalone
            #--template=$(TEMPLATE) \

# Generate all html if the org file change or the template change
$(DST_DIR)/%.html: $(SRC_DIR)/%.md $(TEMPLATE); mkdir -p $(dir $@); $(PANDOC) $< --output $@

# make deploy will deploy the files to my website write your own script
deploy: $(ALL); engine/deploy.sh

website: $(ALL)

.PHONY: clean

clean: ; -rm -rf $(DST_DIR)/*

