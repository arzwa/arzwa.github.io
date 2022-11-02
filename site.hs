{-# LANGUAGE OverloadedStrings #-}
import Data.Monoid (mappend)
import Hakyll
import Text.Pandoc.Options

main :: IO ()
main = hakyll $ do
    match "pdfs/*" $ do
        route   idRoute
        compile copyFileCompiler
    
    match "img/*" $ do
        route   idRoute
        compile copyFileCompiler
    
    match "img/*/*" $ do
        route   idRoute
        compile copyFileCompiler

    match "css/*" $ do
        route   idRoute
        compile compressCssCompiler

    match "posts/*/img/*" $ do
        route   idRoute
        compile copyFileCompiler
    
    match "posts/*/*.svg" $ do
        route   idRoute
        compile copyFileCompiler
    
    match "drafts/*/*.svg" $ do
        route   idRoute
        compile copyFileCompiler
    
    match "drafts/*/img/*" $ do
        route   idRoute
        compile copyFileCompiler

    match (fromList ["about.md", "other.md"]) $ do
        route   $ setExtension "html"
        compile $ pandocMathCompiler
            >>= loadAndApplyTemplate "templates/default.html" defaultContext
            >>= relativizeUrls

    -- should get bibiography in there...
    match mdPostPattern $ do
        route $ setExtension "html"
        compile $ pandocMathCompiler
            >>= loadAndApplyTemplate "templates/post.html"    postCtx
            >>= loadAndApplyTemplate "templates/default.html" postCtx
            >>= relativizeUrls
    
    match mdDraftPattern $ do
        route $ setExtension "html"
        compile $ pandocMathCompiler
            >>= loadAndApplyTemplate "templates/draft.html"    draftCtx
            >>= loadAndApplyTemplate "templates/default.html" draftCtx
            >>= relativizeUrls
    
    create ["archive.html"] $ do
        route idRoute
        compile $ do
            posts <- recentFirst =<< loadAll mdPostPattern
            let archiveCtx =
                    listField "posts" postCtx (return posts) `mappend`
                    constField "title" ""                    `mappend`
                    defaultContext
            makeItem ""
                >>= loadAndApplyTemplate "templates/archive.html" archiveCtx
                >>= loadAndApplyTemplate "templates/default.html" archiveCtx
                >>= relativizeUrls

    match "index.html" $ do
        route idRoute
        compile $ do
            posts <- recentFirst =<< loadAll mdPostPattern
            let indexCtx =
                    listField "posts" postCtx (return posts) `mappend`
                    defaultContext

            getResourceBody
                >>= applyAsTemplate indexCtx
                >>= loadAndApplyTemplate "templates/default.html" indexCtx
                >>= relativizeUrls

    match "templates/*" $ compile templateBodyCompiler

mdPostPattern :: Pattern
mdPostPattern = "posts/*.md" .||. "posts/*/*.md"

mdDraftPattern :: Pattern
mdDraftPattern = "drafts/*.md" .||. "drafts/*/*.md"

postCtx :: Context String
postCtx =
    dateField "date" "%B %e, %Y" `mappend`
    defaultContext

draftCtx :: Context String
draftCtx  =
    dateField "date" "..." `mappend`
    defaultContext

pandocMathCompiler =
    let readerOptions = defaultHakyllReaderOptions
        writerOptions = defaultHakyllWriterOptions
                        { writerHTMLMathMethod = MathJax "" }
    in pandocCompilerWith readerOptions writerOptions
