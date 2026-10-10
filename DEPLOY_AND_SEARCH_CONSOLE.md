# AntennaSim SEO deployment pack

This pack contains **only** SEO-related website files; it does not contain or change your `em-simulator.js`. The replacement `index.html` starts from the October 10 `index(2).html` copy you shared, keeping the original script path (`em-simulator.js?v=20261010-4`), p5.js dependency, startup error handling, and Cloudflare Analytics token intact. If the repository's live `index.html` has changed since that copy, merge the changes into the live file instead of overwriting it.

## Files

- `index.html` — metadata, Schema.org WebApplication JSON-LD, crawlable introduction, small footer About disclosure.
- `about.html` — standalone indexable content and simulator guide; works without JavaScript.
- `robots.txt` — permits crawling and announces the sitemap. Merge with existing robots rules if any.
- `sitemap.xml` — contains the homepage and About page.

## Deploy

1. Commit/upload all four files to the root of the GitHub Pages branch for `www.antennasim.com`.
2. Keep the existing `em-simulator.js` file in the same directory. Do not replace it with an old copy.
3. Confirm the site loads, all tools/presets still work, About opens in the footer, and Cloudflare continues receiving visits.
4. Open `https://www.antennasim.com/about.html`, `/sitemap.xml`, and `/robots.txt` in a browser; they should return the expected HTML/XML/text rather than an error.
5. Confirm that `https://antennasim.com/` redirects to `https://www.antennasim.com/`, or else use one consistent canonical URL and redirect the alternative. This was not verified automatically.

## Google Search Console

1. Visit https://search.google.com/search-console/ and choose **Add property** → **Domain**.
2. Enter `antennasim.com` (no `https://` and no `www`).
3. Copy the unique DNS TXT verification value supplied by Google. Add it to the **authoritative DNS provider** for your domain (GoDaddy if it still hosts DNS; otherwise the provider shown by your nameservers). The host is commonly `@`.
4. Return to Search Console and press **Verify** after DNS propagation.
5. Open **Sitemaps** and submit `https://www.antennasim.com/sitemap.xml`.
6. Use **URL inspection** for `https://www.antennasim.com/` and `https://www.antennasim.com/about.html`. Inspect whether each is indexed, run **Test live URL**, and use **Request indexing** if needed.
7. Recheck the **Pages** and **Performance** reports after Google has crawled the changes. Indexing and ranking are not guaranteed.

The unique Google TXT value cannot be generated here; Search Console issues it to the verified account. I cannot claim the domain is registered or verified until those steps have been completed.
