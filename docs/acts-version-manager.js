/**
 * ACTS Version Manager
 *
 * Handles dynamic version selection and PR build detection for ACTS documentation.
 *
 * Features:
 * - Detects PR builds via pr.json file
 * - Displays PR banner for PR builds
 * - Loads external version selector from deploy repo
 * - Falls back to minimal selector if external JS unavailable
 */

class ActsVersionManager {
  constructor() {
    this.siteRoot = null; // Will be set during discovery
    this.prMetadata = null;
    console.log('[ACTS] Version manager initialized');
  }

  /**
   * Calculate relative path to site root based on current URL depth
   * @returns {string} Relative path (e.g., '.', '..', '../..')
   */
  calculateSiteRoot() {
    const pathname = window.location.pathname;

    // Remove leading/trailing slashes
    const cleanPath = pathname.replace(/^\/+|\/+$/g, '');

    if (!cleanPath) {
      return '.'; // Root path
    }

    // Count directory levels (excluding filename)
    const parts = cleanPath.split('/');
    const depth = parts.length - 1; // Subtract 1 for filename

    if (depth === 0) {
      return '.'; // At root (index.html)
    } else {
      return '../'.repeat(depth).slice(0, -1); // Remove trailing slash
    }
  }

  /**
   * Discover the actual docs root by trying to find index.html at different levels
   * This handles cases where docs are deployed at nested paths (e.g., PR previews)
   * @returns {Promise<string>} Path to docs root
   */
  async discoverDocsRoot() {
    const pathname = window.location.pathname;

    // Special case: if we're already on index.html, we're at the root
    if (pathname.endsWith('/index.html') || pathname.endsWith('/')) {
      console.log('[ACTS] Already at docs root (index.html)');
      return '.';
    }

    const cleanPath = pathname.replace(/^\/+|\/+$/g, '');
    const parts = cleanPath.split('/');
    const maxDepth = parts.length - 1; // Don't include the filename

    // Try each level from current directory up to domain root
    for (let depth = 0; depth <= maxDepth; depth++) {
      const path = depth === 0 ? '.' : '../'.repeat(depth).slice(0, -1);

      try {
        const controller = new AbortController();
        const timeout = setTimeout(() => controller.abort(), 500);

        const response = await fetch(`${path}/index.html`, {
          method: 'HEAD', // Just check existence, don't download
          signal: controller.signal
        });

        clearTimeout(timeout);

        if (response.ok) {
          console.log(`[ACTS] Discovered docs root at: ${path}`);
          return path;
        }
      } catch (error) {
        // Continue searching
      }
    }

    // Fallback to calculated site root if discovery fails
    console.warn('[ACTS] Could not discover docs root, using calculated path');
    return this.calculateSiteRoot();
  }

  /**
   * Check for PR build by attempting to fetch pr.json
   * @returns {Promise<object|null>} PR metadata or null
   */
  async checkForPRBuild() {
    try {
      const controller = new AbortController();
      const timeout = setTimeout(() => controller.abort(), 2000);

      const response = await fetch(`${this.siteRoot}/pr.json`, {
        cache: 'no-store',
        signal: controller.signal
      });

      clearTimeout(timeout);

      if (!response.ok) {
        console.log('[ACTS] No pr.json found - not a PR build');
        return null;
      }

      const data = await response.json();

      // Validate structure
      if (!data.pr || !data.url) {
        console.warn('[ACTS] Invalid pr.json structure - missing required fields');
        return null;
      }

      console.log('[ACTS] PR build detected:', data);
      return data;
    } catch (error) {
      if (error.name === 'AbortError') {
        console.warn('[ACTS] pr.json fetch timeout');
      } else {
        console.log('[ACTS] No pr.json found (expected for non-PR builds)');
      }
      return null;
    }
  }

  /**
   * Render PR build banner at the top of #doc-content
   * @param {object} metadata - PR metadata from pr.json
   */
  renderPRBanner(metadata) {
    const docContent = document.querySelector('#doc-content');
    if (!docContent) {
      console.warn('[ACTS] Cannot render PR banner - #doc-content not found');
      return;
    }

    const banner = document.createElement('div');
    banner.id = 'acts-pr-banner';

    // Build banner content
    let content = `<strong>PR Build #${metadata.pr}</strong>`;

    if (metadata.branch) {
      content += ` from branch <code>${metadata.branch}</code>`;
    }

    if (metadata.title) {
      content += ` - ${metadata.title}`;
    }

    if (metadata.url) {
      content = `${content} (<a href="${metadata.url}" target="_blank" rel="noopener">View PR</a>)`;
    }

    banner.innerHTML = content;

    // Insert at the top of #doc-content
    docContent.prepend(banner);

    console.log('[ACTS] PR banner rendered at top of #doc-content');
  }

  /**
   * Load external version-selector.js from site root
   * @returns {Promise<boolean>} True if loaded successfully
   */
  async loadExternalVersionSelector() {
    return new Promise((resolve) => {
      const script = document.createElement('script');
      script.src = `${this.siteRoot}/version-selector.js`;

      const timeout = setTimeout(() => {
        console.warn('[ACTS] External version-selector.js timeout after 5s');
        script.remove();
        this.renderMinimalSelector();
        resolve(false);
      }, 5000);

      script.onload = () => {
        clearTimeout(timeout);

        if (typeof window.initVersionSelector === 'function') {
          try {
            console.log('[ACTS] External version-selector.js loaded, initializing...');
            window.initVersionSelector({
              container: '#acts-version-selector',
              currentPath: window.location.pathname,
              siteRoot: this.siteRoot
            });
            console.log('[ACTS] Version selector initialized successfully');
            resolve(true);
          } catch (error) {
            console.error('[ACTS] version-selector.js init failed:', error);
            this.renderMinimalSelector();
            resolve(false);
          }
        } else {
          console.warn('[ACTS] version-selector.js missing initVersionSelector function');
          this.renderMinimalSelector();
          resolve(false);
        }
      };

      script.onerror = () => {
        clearTimeout(timeout);
        console.log('[ACTS] External version-selector.js not available (expected for local builds)');
        this.renderMinimalSelector();
        resolve(false);
      };

      document.head.appendChild(script);
    });
  }

  /**
   * Render minimal fallback version selector
   */
  renderMinimalSelector() {
    // Disabled for now - remove the empty container instead
    const container = document.querySelector('#acts-version-selector');
    if (container) {
      container.remove();
      console.log('[ACTS] Minimal selector disabled - removed empty container');
    }
  }

  /**
   * Main initialization
   */
  async init() {
    // Discover the actual docs root (handles nested deployments)
    this.siteRoot = await this.discoverDocsRoot();
    console.log('[ACTS] Using docs root:', this.siteRoot);

    // Check for PR build
    this.prMetadata = await this.checkForPRBuild();

    if (this.prMetadata) {
      // This is a PR build - show banner only
      this.renderPRBanner(this.prMetadata);
    } else {
      // Regular build - show version selector
      // Create container first to avoid layout shift
      const sideNav = document.querySelector('#side-nav');
      if (sideNav) {
        const container = document.createElement('div');
        container.id = 'acts-version-selector';
        sideNav.prepend(container);

        // Load external version selector
        await this.loadExternalVersionSelector();
      } else {
        console.warn('[ACTS] Cannot create version selector - #side-nav not found');
      }
    }
  }
}

// Initialize when DOM ready (jQuery is already loaded by Doxygen)
$(document).ready(function() {
  const versionManager = new ActsVersionManager();
  versionManager.init();
});
