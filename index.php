<?php
$search = $_GET['search'] ?? '';
$cwd = getcwd();
$relPath = $_SERVER['REQUEST_URI'];
$pathParts = array_filter(explode('/', $relPath));

$subdirs = array_filter(glob('*'), 'is_dir');
$files = array_filter(glob('*'), 'is_file');

function is_plot($f) {
    return preg_match('/\.(png|jpe?g|gif|pdf)$/i', $f);
}
function matches_search($f, $search) {
    return empty($search) || stripos($f, $search) !== false;
}


$entries = glob('*', GLOB_NOSORT);

// Directories sorted by mtime (desc)
$subdirs = array_values(array_filter($entries, 'is_dir'));
usort($subdirs, function($a, $b) {
    $mb = @filemtime($b) ?: 0;
    $ma = @filemtime($a) ?: 0;
    if ($mb === $ma) return strcasecmp($a, $b);  // tie-break by name
    return $mb <=> $ma;                           // newest first
});

// Files (keep if you also want sorted files; otherwise remove this part)
$files = array_values(array_filter($entries, 'is_file'));
usort($files, function($a, $b) {
    $mb = @filemtime($b) ?: 0;
    $ma = @filemtime($a) ?: 0;
    if ($mb === $ma) return strcasecmp($a, $b);
    return $mb <=> $ma;
});
?>
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title><?php echo basename($cwd); ?></title>
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/css/bootstrap.min.css" rel="stylesheet">
  <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap-icons@1.11.1/font/bootstrap-icons.css">
  <style>
    body > .container-fluid { margin-top: 20px; padding: 0 20px; }
    .empty-text { font-style: italic; font-size: 0.9rem; }
	/* Card size is driven by a CSS variable so the JS slider resizes the whole box */
	#plot-listing .card {
	  margin: 5px;
	  width: var(--plot-card-width, 240px);
	  max-width: calc(100vw - 60px); /* never wider than the page: 1 plot per page at max */
	}
	.card-img-top {
	  width: 100%;
	  height: auto;
	  object-fit: contain;
	}
	/* JSROOT histogram drawing area scales with the card (4:3 aspect ratio) */
	.jsroot-plot {
	  width: 100%;
	  aspect-ratio: 4 / 3;
	  background: #fff;
	}
  </style>
</head>
<body>
  <!-- Breadcrumb navigation -->
  <nav class="navbar navbar-light bg-light border-bottom">
    <div class="container-fluid d-flex justify-content-between align-items-center">
	<div class="d-flex align-items-center">
	  <a href="/" class="text-decoration-none" style="color: blue;"><i class="bi bi-house-door"></i></a>

	  <?php
	    $uri = parse_url($_SERVER['REQUEST_URI'], PHP_URL_PATH);
	    $parts = explode('/', trim($uri, '/'));
	    $accum = '';
	    foreach ($parts as $i => $part) {
		$accum .= '/' . $part;
		echo '<span class="mx-2">/</span>';
		if ($i < count($parts) - 1) {
		    echo '<a href="' . htmlspecialchars($accum) . '" class="text-decoration-none" style="color: blue;">' . htmlspecialchars($part) . '</a>';
		} else {
		    echo '<span class="fw-light" style="color: blue;">' . htmlspecialchars($part) . '</span>';
		}
	    }
	  ?>
	</div>
      <form class="d-flex" method="get">
        <input class="form-control me-2" type="search" name="search" placeholder="Pattern(s)" value="<?php echo htmlspecialchars($search); ?>">
        <button class="btn btn-outline-success" type="submit">Search</button>
      </form>
    </div>
  </nav>

  <div class="container-fluid">
    <!-- Subdirectories -->
    <h5 class="mt-4">Directories</h5>
    <?php if (!empty($subdirs)): ?>
      <ul>
        <?php foreach ($subdirs as $dir): ?>
          <?php if (matches_search($dir, $search)): ?>
            <li><a href="<?php echo htmlspecialchars($dir); ?>"><?php echo htmlspecialchars($dir); ?></a></li>
          <?php endif; ?>
        <?php endforeach; ?>
      </ul>
    <?php else: ?>
      <p class="empty-text">No directories found.</p>
    <?php endif; ?>

	<!-- Plots -->
	<h5 class="mt-4">Plots</h5>

	<!-- Image size slider -->
	<div class="d-flex align-items-center mb-2" style="gap: 10px;">
  		<label for="imgSize" class="mb-0">Image size:</label>
  		<input
    		type="range"
    		id="imgSize"
    		min="120"
    		max="1600"
    		value="240"
    		style="width: 220px;"
  		>
		</div>

	<div class="d-flex flex-wrap" id="plot-listing">
	<?php
	$displayed = [];
	$formats = ['png', 'pdf', 'root', 'C', 'jpg', 'jpeg', 'gif']; // include C and others

	foreach ($files as $file) {
	    if (!is_plot($file) && !preg_match('/\.(C)$/i', $file)) continue;
	    if (!matches_search($file, $search)) continue;

	    $ext = pathinfo($file, PATHINFO_EXTENSION);
	    $base = preg_replace('/\.(png|jpe?g|gif|pdf|root|C)$/i', '', $file);

	    // One box per plot: names that differ only by a leading "c_" are the same plot
	    $canon = preg_replace('/^c_/', '', $base);

	    if (in_array($canon, $displayed)) continue;



      $available = [];

      // Formati visualizzabili/scaricabili
      foreach (['png', 'pdf', 'C', 'jpg', 'jpeg', 'gif'] as $fmt) {
          foreach ([$canon, 'c_' . $canon] as $b) {
              $candidate = $b . '.' . $fmt;

              if (isset($available[$fmt])) {
                  continue;
              }

              if (file_exists($candidate)) {
                  $available[$fmt] = $candidate;
              }
          }
      }

      // Se non c'è nessun formato visualizzabile, non mostrare nulla.
      // In particolare: un .root da solo NON genera una card.
      if (empty($available)) {
          continue;
      }

      // Aggiungi ROOT solo se esiste e solo se esiste già almeno
      // un formato visualizzabile.
      foreach ([$canon, 'c_' . $canon] as $b) {
          $candidate = $b . '.root';

          if (file_exists($candidate)) {
              $available['root'] = $candidate;
              break;
          }
      }


	    // Show only if there's an image to display
	    $thumb = $available['png'] ?? null;
	    if (!$thumb) continue;

	    $displayed[] = $canon;
	    $baseName = basename($canon);

	    echo '<div class="card">';
	    echo '<a href="' . htmlspecialchars($thumb) . '" target="_blank">';
	    echo '<div class="card-header text-danger fw-bold text-center">' . htmlspecialchars($baseName) . '</div>';
	    echo '<img src="' . htmlspecialchars($thumb) . '" class="card-img-top" alt="' . $thumb . '">';
	    echo '</a>';

	    echo '<div class="card-footer text-center">';
	    foreach ($available as $fmt => $path) {
		$color = match(strtolower($fmt)) {
		    'png' => 'primary',
		    'pdf' => 'danger',
		    'root' => 'dark',
		    'c' => 'success',
		    'jpg', 'jpeg' => 'warning',
		    'gif' => 'info',
		    default => 'secondary',
		};
		echo '<a href="' . htmlspecialchars($path) . '" target="_blank" class="badge bg-' . $color . ' mx-1 text-decoration-none">' . htmlspecialchars($fmt) . '</a>';
	    }
	    echo '</div></div>';
	}
	if (empty($displayed)) {
	    echo "<p class='empty-text' id='no-plots'>No plots to display</p>";
	}

	// .root files (matching search) whose histograms will be drawn client-side with JSROOT
	$rootFiles = array_values(array_filter($files, function($f) use ($search) {
	    return preg_match('/\.root$/i', $f) && matches_search($f, $search);
	}));
	// Canonical names of plots already shown as images, so JSROOT won't duplicate them
	$displayedCanon = array_values(array_map('basename', $displayed));
	?>
	</div>

    <!-- Other Files -->
    <h5 class="mt-4">Other files</h5>
    <ul>
      <?php
      // List every non-image file (e.g. .root, .C, .txt, ...) except index.php
      $others = array_filter($files, function($f) use ($search) {
          return !is_plot($f) && matches_search($f, $search) && basename($f) !== 'index.php';
      });
      if (!empty($others)):
        foreach ($others as $f):
            echo "<li><a href=\"" . htmlspecialchars($f) . "\">" . htmlspecialchars($f) . "</a></li>";
        endforeach;
      else:
        echo "<p class='empty-text'>No files to display</p>";
      endif;
      ?>
    </ul>

    <!-- Scroll to Top -->
    <div class="text-start mt-3">
      <button id="scroll-top" class="btn btn-outline-primary btn-sm">To top</button>
    </div>
  </div>
	<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/js/bootstrap.bundle.min.js"></script>
	<script>
	  document.getElementById('scroll-top').onclick =
	    () => window.scrollTo({top: 0, behavior: 'smooth'});
	
	  // Image size slider logic: resizes the whole plot box (card)
	  const imgSizeSlider = document.getElementById('imgSize');
	  function updateImageSize() {
	    const w = imgSizeSlider.value;
	    document.documentElement.style.setProperty('--plot-card-width', w + 'px');
	  }
	  if (imgSizeSlider) {
	    imgSizeSlider.addEventListener('input', updateImageSize);
	    updateImageSize(); // initialise on load
	  }
	</script>

</body>
</html>
