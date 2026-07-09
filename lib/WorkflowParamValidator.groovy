import java.util.regex.Pattern

class WorkflowParamValidator {
    private static final Pattern INTEGER = Pattern.compile(/[0-9]+/)
    private static final Pattern TOKEN = Pattern.compile(/[A-Za-z0-9][A-Za-z0-9._+-]*/)
    private static final Pattern PATH_VALUE = Pattern.compile(/[A-Za-z0-9._+,:=@%\/-]+/)
    private static final Pattern HEADER_VALUE = Pattern.compile(/[^\p{Cntrl}$`"';|&<>]+/)
    private static final Set DROPLET_PROTOCOLS = ['10xv1', '10xv1a', '10xv1i', '10xv2', '10xv3', 'drop-seq', 'seq-well', '10x5prime'] as Set
    private static final Set SMART_PROTOCOLS = ['smart-seq', 'smart-seq2', 'smarter', 'smart-like'] as Set
    private static final Set TERTIARY_WORKFLOWS = ['none', 'scanpy-workflow', 'scanpy-NFworkflow'] as Set
    private static final Set YES_NO = ['yes', 'no'] as Set

    static void validate(def params) {
        requirePath(params, 'resultsRoot')
        requirePath(params, 'masterWorkflow')
        requireToken(params, 'expName')
        requireProtocolList(params)

        [
            'rawMatrix',
            'referenceFasta',
            'referenceGtf',
            'geneMetadata',
            'cellMetadata',
            'condensedSdrf',
            'projectFile'
        ].each { requirePath(params, it) }

        requireInteger(params, 'largeMatrixThreshold')
        requireInteger(params, 'topmarkersForSummary')
        requireFields(params.fields)

        if (has(params, 'tertiary')) {
            requireEnum(params, 'tertiary', YES_NO)
        }

        def tertiaryWorkflow = has(params, 'tertiaryWorkflow') ? params.tertiaryWorkflow.toString() : 'none'
        if (!(tertiaryWorkflow in TERTIARY_WORKFLOWS)) {
            throw new IllegalArgumentException("params.tertiaryWorkflow must be one of ${TERTIARY_WORKFLOWS}; got '${tertiaryWorkflow}'")
        }
        if (tertiaryWorkflow in ['scanpy-workflow', 'scanpy-NFworkflow']) {
            [
                'rawFilteredMatrix',
                'normalisedMatrix',
                'clusters',
                'tsneDir',
                'umapDir',
                'markersDir'
            ].each { requirePath(params, it) }
        }
        if (tertiaryWorkflow == 'scanpy-NFworkflow') {
            requirePath(params, 'tertiarySoftwareReport')
        } else {
            optionalPath(params, 'tertiarySoftwareReport')
        }

        optionalPath(params, 'tpmMatrix')
    }

    static String shellQuote(value) {
        "'" + value.toString().replace("'", "'\"'\"'") + "'"
    }

    private static void requireProtocolList(def params) {
        requireValue(params, 'params.protocolList', 'protocolList')
        def protocols = params.protocolList.toString().split(',') as List
        if (protocols.isEmpty()) {
            throw new IllegalArgumentException('params.protocolList must contain at least one protocol')
        }
        protocols.each { protocol ->
            if (!((protocol in DROPLET_PROTOCOLS) || (protocol in SMART_PROTOCOLS))) {
                throw new IllegalArgumentException("Unsupported protocol in params.protocolList: '${protocol}'")
            }
        }
    }

    private static void requireFields(def fields) {
        if (fields == null) {
            throw new IllegalArgumentException('Missing workflow params.fields settings')
        }
        requireNestedHeader(fields, 'params.fields', 'run')
        optionalNestedHeader(fields, 'params.fields', 'techrep')
    }

    private static void requirePath(def params, String name) {
        requireValue(params, "params.${name}", name)
        assertPattern(params.get(name), "params.${name}", PATH_VALUE)
    }

    private static void optionalPath(def params, String name) {
        if (has(params, name) && params.get(name) != null && params.get(name).toString() != '') {
            assertPattern(params.get(name), "params.${name}", PATH_VALUE)
        }
    }

    private static void requireToken(def params, String name) {
        requireValue(params, "params.${name}", name)
        assertPattern(params.get(name), "params.${name}", TOKEN)
    }

    private static void requireInteger(def params, String name) {
        requireValue(params, "params.${name}", name)
        assertPattern(params.get(name), "params.${name}", INTEGER)
    }

    private static void requireEnum(def params, String name, Set allowed) {
        requireValue(params, "params.${name}", name)
        def text = params.get(name).toString()
        if (!(text in allowed)) {
            throw new IllegalArgumentException("params.${name} must be one of ${allowed}; got '${text}'")
        }
    }

    private static void requireNestedHeader(def params, String scope, String name) {
        requireValue(params, "${scope}.${name}", name)
        assertPattern(params.get(name), "${scope}.${name}", HEADER_VALUE)
    }

    private static void optionalNestedHeader(def params, String scope, String name) {
        if (has(params, name) && params.get(name) != null && params.get(name).toString() != '') {
            assertPattern(params.get(name), "${scope}.${name}", HEADER_VALUE)
        }
    }

    private static void requireValue(def params, String label, String key) {
        if (params == null || !has(params, key) || params.get(key) == null || params.get(key).toString() == '') {
            throw new IllegalArgumentException("Missing required workflow parameter ${label}")
        }
    }

    private static void assertPattern(value, String label, Pattern pattern) {
        def text = value == null ? '' : value.toString()
        if (!pattern.matcher(text).matches()) {
            throw new IllegalArgumentException("Invalid workflow parameter ${label}: '${text}'")
        }
    }

    private static boolean has(def params, String key) {
        params != null && params.containsKey(key)
    }
}
